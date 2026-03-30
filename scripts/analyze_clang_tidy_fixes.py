import yaml, collections

data = yaml.safe_load(open("fixes.yaml"))
counts = collections.Counter(d["DiagnosticName"] for d in data.get("Diagnostics", []))
for name, n in counts.most_common():
    print(f"{n:4d}  {name}")
