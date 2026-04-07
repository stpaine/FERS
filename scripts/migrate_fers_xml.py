#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import argparse
import sys


def indent(elem, level=0):
    """Fallback indentation for Python versions < 3.9"""
    i = "\n" + level * "    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "    "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for child in elem:
            indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def get_user_input(prompt_text, default_val):
    """Gets user input with a fallback for non-interactive environments."""
    try:
        user_input = input(prompt_text).strip()
        return user_input if user_input else default_val
    except EOFError:
        print(default_val)
        return default_val


def resolve_duplicate_names(root):
    """
    Scans the XML tree for duplicate names, prompts the user to rename them,
    and updates any references to those names to maintain original usage.
    """
    definitions = {}
    references = {}

    # Pass 1: Collect all named elements and all references
    for el in root.iter():
        name = el.get("name")
        if name is not None:
            definitions.setdefault(name, []).append(el)

        # Track attributes that reference other elements
        for attr in ["pulse", "timing", "antenna"]:
            val = el.get(attr)
            if val is not None:
                references.setdefault(val, []).append((el, attr))

    used_names = set(definitions.keys())

    # Pass 2: Resolve duplicates
    for old_name, elements in definitions.items():
        if len(elements) > 1:
            print(f"\n[!] Duplicate name detected: '{old_name}' is used by {len(elements)} elements.")
            used_names.remove(old_name)  # Will be replaced entirely
            new_names = []

            for i, el in enumerate(elements):
                tag = el.tag
                counter = i + 1
                default_name = f"{old_name}{counter}"

                # Ensure default name isn't already taken
                while default_name in used_names:
                    counter += 1
                    default_name = f"{old_name}{counter}"

                prompt = f"    Rename <{tag}> (occurrence {i + 1}) [{default_name}]: "
                new_name = get_user_input(prompt, default_name)

                # Ensure user's chosen name isn't already taken
                while new_name in used_names:
                    print(f"    Error: Name '{new_name}' is already in use. Choose another.")
                    new_name = get_user_input(prompt, default_name)

                el.set("name", new_name)
                used_names.add(new_name)
                new_names.append(new_name)

            # Pass 3: Update references to the old name
            refs = references.get(old_name, [])
            if refs:
                print(f"\n    [!] There are {len(refs)} references to '{old_name}'.")
                for ref_el, attr in refs:
                    ref_tag = ref_el.tag
                    ref_name = ref_el.get("name", "<unnamed>")
                    default_ref = new_names[0]

                    print(f"    Element <{ref_tag} name='{ref_name}'> references {attr}='{old_name}'.")
                    print(f"    Available options: {', '.join(new_names)}")

                    prompt = f"    Update reference to [{default_ref}]: "
                    chosen_ref = get_user_input(prompt, default_ref)

                    while chosen_ref not in new_names:
                        print(f"    Invalid choice. Must be one of: {', '.join(new_names)}")
                        chosen_ref = get_user_input(prompt, default_ref)

                    ref_el.set(attr, chosen_ref)


def extract_and_append(old_parent, new_parent, old_tag, new_tag=None):
    """Finds an element in the old parent and appends it to the new parent.
       This enforces strict ordering by pulling elements exactly when called."""
    if new_tag is None:
        new_tag = old_tag
    el = old_parent.find(old_tag)
    if el is not None:
        new_el = ET.SubElement(new_parent, new_tag)
        new_el.text = el.text
        for k, v in el.attrib.items():
            new_el.set(k, v)
        return new_el
    return None


def copy_attributes(old_el, new_el, rename_map=None, drop_list=None):
    """Copies attributes from old element to new element with optional renaming/dropping."""
    if rename_map is None: rename_map = {}
    if drop_list is None: drop_list = []
    for k, v in old_el.attrib.items():
        if k in drop_list:
            continue
        new_k = rename_map.get(k, k)
        new_el.set(new_k, v)


def migrate_xml(input_file, output_file):
    try:
        tree = ET.parse(input_file)
        old_root = tree.getroot()
    except ET.ParseError as e:
        print(f"Error parsing XML file: {e}")
        sys.exit(1)

    if old_root.tag != "simulation":
        print("Error: Root element is not <simulation>.")
        sys.exit(1)

    # Pre-process the XML to resolve duplicate names interactively
    resolve_duplicate_names(old_root)

    new_root = ET.Element("simulation")
    new_root.set("name", old_root.get("name", "migrated_simulation"))

    # 1. Process <parameters> (Must be first)
    old_params = old_root.find("parameters")
    if old_params is not None:
        new_params = ET.SubElement(new_root, "parameters")
        # Enforce strict XSD sequence order
        extract_and_append(old_params, new_params, "starttime")
        extract_and_append(old_params, new_params, "endtime")
        extract_and_append(old_params, new_params, "rate")
        extract_and_append(old_params, new_params, "c")
        extract_and_append(old_params, new_params, "interprate", "simSamplingRate")
        extract_and_append(old_params, new_params, "randomseed")
        extract_and_append(old_params, new_params, "adc_bits")
        extract_and_append(old_params, new_params, "oversample")

        if old_params.find("export") is not None:
            print("\nNotice: <export> element is deprecated and was removed.")

    # 2. Process other root children
    for child in old_root:
        if child.tag == "parameters":
            continue

        # <pulse> -> <waveform>
        elif child.tag == "pulse":
            new_wave = ET.SubElement(new_root, "waveform")
            new_wave.set("name", child.get("name"))

            # Enforce strict XSD sequence order
            extract_and_append(child, new_wave, "power")
            extract_and_append(child, new_wave, "carrier", "carrier_frequency")

            ptype = child.get("type")
            if ptype == "continuous" or ptype == "cw":
                print(
                    f"\nWarning: <pulse name='{child.get('name')}'> has type='continuous'. FERS now has a native CW mode.")
                print("\t Please consult the documentation on simulating continuous mode.")

            if ptype == "file":
                pf = ET.SubElement(new_wave, "pulsed_from_file")
                pf.set("filename", child.get("filename", ""))
            else:
                ET.SubElement(new_wave, "cw")

        # <timing>
        elif child.tag == "timing":
            new_timing = ET.SubElement(new_root, "timing")
            new_timing.set("name", child.get("name"))

            # synconpulse default changed from true to false. Preserve old default.
            synconpulse = child.get("synconpulse", "true")
            new_timing.set("synconpulse", synconpulse)

            # Enforce strict XSD sequence order
            extract_and_append(child, new_timing, "frequency")
            extract_and_append(child, new_timing, "freq_offset")
            extract_and_append(child, new_timing, "random_freq_offset", "random_freq_offset_stdev")
            extract_and_append(child, new_timing, "phase_offset")
            extract_and_append(child, new_timing, "random_phase_offset", "random_phase_offset_stdev")

            for ne in child.findall("noise_entry"):
                new_ne = ET.SubElement(new_timing, "noise_entry")
                extract_and_append(ne, new_ne, "alpha")
                extract_and_append(ne, new_ne, "weight")

        # <antenna>
        elif child.tag == "antenna":
            new_ant = ET.SubElement(new_root, "antenna")
            new_ant.set("name", child.get("name"))
            new_ant.set("pattern", child.get("pattern"))
            if child.get("filename"):
                new_ant.set("filename", child.get("filename"))

            if child.get("module") or child.get("function"):
                print(
                    f"\nWarning: Python antennas (module/function) are deprecated. Removed from antenna '{child.get('name')}'.")

            # Enforce strict XSD sequence order
            for prop in ["alpha", "beta", "gamma", "diameter", "azscale", "elscale", "efficiency"]:
                extract_and_append(child, new_ant, prop)

        # <platform>
        elif child.tag == "platform":
            new_plat = ET.SubElement(new_root, "platform")
            new_plat.set("name", child.get("name"))

            # 2.1 Motion path (Must be first in platform)
            old_mp = child.find("motionpath")
            if old_mp is not None:
                new_mp = ET.SubElement(new_plat, "motionpath")
                interp = old_mp.get("interpolation", "static")
                if interp == "python":
                    print(
                        f"\nWarning: Python interpolation is deprecated. Defaulting to 'static' for platform '{child.get('name')}'.")
                    interp = "static"
                new_mp.set("interpolation", interp)

                for wp in old_mp.findall("positionwaypoint"):
                    new_wp = ET.SubElement(new_mp, "positionwaypoint")
                    # Enforce strict XSD sequence order
                    extract_and_append(wp, new_wp, "x")
                    extract_and_append(wp, new_wp, "y")
                    extract_and_append(wp, new_wp, "altitude")
                    extract_and_append(wp, new_wp, "time")

            # 2.2 Rotation paths (Must be second in platform)
            old_rp = child.find("rotationpath")
            old_fr = child.find("fixedrotation")

            if old_rp is not None:
                new_rp = ET.SubElement(new_plat, "rotationpath")
                copy_attributes(old_rp, new_rp)
                for rw in old_rp.findall("rotationwaypoint"):
                    new_rw = ET.SubElement(new_rp, "rotationwaypoint")
                    # Enforce strict XSD sequence order
                    extract_and_append(rw, new_rw, "azimuth")
                    extract_and_append(rw, new_rw, "elevation")
                    extract_and_append(rw, new_rw, "time")
            elif old_fr is not None:
                new_fr = ET.SubElement(new_plat, "fixedrotation")
                # Enforce strict XSD sequence order
                extract_and_append(old_fr, new_fr, "startazimuth")
                extract_and_append(old_fr, new_fr, "startelevation")
                extract_and_append(old_fr, new_fr, "azimuthrate")
                extract_and_append(old_fr, new_fr, "elevationrate")

            # 2.3 Radar nodes (Must be after paths)
            for rn in child:
                if rn.tag in ["monostatic", "transmitter", "receiver"]:
                    new_rn = ET.SubElement(new_plat, rn.tag)

                    if rn.get("type") == "continuous" or rn.get("type") == "cw":
                        print(
                            f"\nWarning: <{rn.tag} name='{rn.get('name')}'> has type='continuous'. FERS now has a native CW mode.")
                        print(
                            "\t Please consult the documentation on simulating continuous mode.")

                    # Receivers do not have a waveform attribute in the new schema
                    if rn.tag == "receiver":
                        copy_attributes(rn, new_rn, drop_list=["type", "pulse"])
                    else:
                        copy_attributes(rn, new_rn, rename_map={"pulse": "waveform"}, drop_list=["type"])

                    is_pulsed = (rn.get("type") == "pulsed" or rn.find("prf") is not None)

                    if is_pulsed:
                        mode_el = ET.SubElement(new_rn, "pulsed_mode")
                        # Enforce strict XSD sequence order
                        extract_and_append(rn, mode_el, "prf")
                        extract_and_append(rn, mode_el, "window_skip")
                        extract_and_append(rn, mode_el, "window_length")
                    else:
                        ET.SubElement(new_rn, "cw_mode")

                    extract_and_append(rn, new_rn, "noise_temp")

                elif rn.tag == "target":
                    new_tgt = ET.SubElement(new_plat, "target")
                    copy_attributes(rn, new_tgt)

                    old_rcs = rn.find("rcs")
                    if old_rcs is not None:
                        new_rcs = ET.SubElement(new_tgt, "rcs")
                        copy_attributes(old_rcs, new_rcs)
                        extract_and_append(old_rcs, new_rcs, "value")

                    old_model = rn.find("model")
                    if old_model is not None:
                        new_model = ET.SubElement(new_tgt, "model")
                        copy_attributes(old_model, new_model)
                        extract_and_append(old_model, new_model, "k")

        # <include>
        elif child.tag == "include":
            new_inc = ET.SubElement(new_root, "include")
            new_inc.text = child.text

        # <multipath>
        elif child.tag == "multipath":
            print("\nNotice: <multipath> element is no longer supported in the new schema and was removed.")

    # Pretty print the XML
    if hasattr(ET, "indent"):
        ET.indent(new_root, space="    ", level=0)
    else:
        indent(new_root)

    # Write to file
    tree = ET.ElementTree(new_root)
    tree.write(output_file, encoding="UTF-8", xml_declaration=True)
    print(f"\nSuccessfully migrated XML to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Migrate FERS XML from the old schema to the new schema.")
    parser.add_argument("input", help="Path to the old XML file")
    parser.add_argument("output", help="Path to save the migrated XML file")

    args = parser.parse_args()
    migrate_xml(args.input, args.output)
