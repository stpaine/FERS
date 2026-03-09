#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "serial/libxml_wrapper.h"

namespace
{
	std::string buildXml(const std::string& content) { return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" + content; }

	std::string writeTempXmlFile(const std::string& content)
	{
		const std::string path = "libxml_wrapper_test.xml";
		std::ofstream out(path, std::ios::binary);
		out << content;
		return path;
	}

	std::vector<unsigned char> toBytes(const std::string& content)
	{
		return std::vector<unsigned char>(content.begin(), content.end());
	}
}

TEST_CASE("XmlElement handles invalid nodes", "[serial][xml]")
{
	XmlElement element(nullptr);
	REQUIRE_FALSE(element.isValid());
	REQUIRE(element.getText().empty());
	REQUIRE_FALSE(element.childElement("child").isValid());
}

TEST_CASE("XmlElement reads and writes text", "[serial][xml]")
{
	XmlDocument doc;
	const XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);
	root.setText("alpha");
	root.setText("beta");
	root.setText("gamma");
	root.setText("delta");

	// Ensure text content is set deterministically and retrieved correctly.
	REQUIRE(root.getText() == "delta");
	REQUIRE(root.name() == "root");
}

TEST_CASE("XmlElement manages attributes", "[serial][xml]")
{
	XmlDocument doc;
	const XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);
	root.setAttribute("mode", "cw");

	REQUIRE(XmlElement::getSafeAttribute(root, "mode") == "cw");
	REQUIRE_THROWS_AS(XmlElement::getSafeAttribute(root, "missing"), XmlException);
}

TEST_CASE("XmlElement handles child indexing and filtering", "[serial][xml]")
{
	XmlDocument doc;
	const XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	const XmlElement alpha0 = root.addChild("alpha");
	alpha0.setText("A0");
	const XmlElement beta0 = root.addChild("beta");
	beta0.setText("B0");
	const XmlElement alpha1 = root.addChild("alpha");
	alpha1.setText("A1");

	REQUIRE(root.childElement("alpha", 0).getText() == "A0");
	REQUIRE(root.childElement("alpha", 1).getText() == "A1");
	REQUIRE(root.childElement("beta", 0).getText() == "B0");
	REQUIRE_FALSE(root.childElement("alpha", 2).isValid());
	REQUIRE_FALSE(root.childElement("gamma", 0).isValid());
}

TEST_CASE("XmlDocument load and dump roundtrip", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><value>42</value></root>");
	REQUIRE(doc.loadString(xml));

	const XmlElement root = doc.getRootElement();
	REQUIRE(root.name() == "root");
	REQUIRE(root.childElement("value", 0).getText() == "42");

	const std::string dumped = doc.dumpToString();
	REQUIRE(dumped.find("<value>42</value>") != std::string::npos);
}

TEST_CASE("XmlDocument loads files and persists output", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><value>17</value></root>");
	const std::string input_path = writeTempXmlFile(xml);

	REQUIRE(doc.loadFile(input_path));
	REQUIRE(doc.getRootElement().childElement("value", 0).getText() == "17");

	const std::string output_path = "libxml_wrapper_test_out.xml";
	REQUIRE(doc.saveFile(output_path));

	std::ifstream in(output_path, std::ios::binary);
	std::string saved((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	REQUIRE(saved.find("<value>17</value>") != std::string::npos);

	std::remove(input_path.c_str());
	std::remove(output_path.c_str());
}

TEST_CASE("XmlDocument requires root element", "[serial][xml]")
{
	XmlDocument doc;
	const XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);
	REQUIRE(doc.getRootElement().name() == "root");
}

TEST_CASE("XmlDocument throws when root is missing", "[serial][xml]")
{
	XmlDocument doc;
	REQUIRE_THROWS_AS(doc.getRootElement(), std::runtime_error);
}

TEST_CASE("XmlDocument throws when setting root on moved-from doc", "[serial][xml]")
{
	XmlDocument doc;
	XmlDocument moved_doc(std::move(doc));
	const XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	moved_doc.setRootElement(root);
	REQUIRE_THROWS_AS(doc.setRootElement(root), std::runtime_error);
}

TEST_CASE("XmlDocument throws when getting root on moved-from doc", "[serial][xml]")
{
	XmlDocument doc;
	XmlDocument moved_doc(std::move(doc));
	const XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	moved_doc.setRootElement(root);
	REQUIRE_THROWS_AS(doc.getRootElement(), std::runtime_error);
}

TEST_CASE("XmlDocument validates with DTD", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><child>7</child></root>");
	REQUIRE(doc.loadString(xml));

	const std::string dtd = "<!ELEMENT root (child)>\n<!ELEMENT child (#PCDATA)>";
	const auto dtd_bytes = toBytes(dtd);
	REQUIRE(doc.validateWithDtd(dtd_bytes));
}

TEST_CASE("XmlDocument fails DTD validation for invalid structure", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><other>7</other></root>");
	REQUIRE(doc.loadString(xml));

	const std::string dtd = "<!ELEMENT root (child)>\n<!ELEMENT child (#PCDATA)>";
	const auto dtd_bytes = toBytes(dtd);
	REQUIRE_THROWS_AS(doc.validateWithDtd(dtd_bytes), XmlException);
}

TEST_CASE("XmlDocument fails with malformed DTD", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><child>7</child></root>");
	REQUIRE(doc.loadString(xml));

	const std::string bad_dtd = "<!ELEMENT root (child";
	const auto dtd_bytes = toBytes(bad_dtd);
	REQUIRE_THROWS_AS(doc.validateWithDtd(dtd_bytes), XmlException);
}

TEST_CASE("XmlDocument validates with XSD", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><child>9</child></root>");
	REQUIRE(doc.loadString(xml));

	const std::string xsd = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
							"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">"
							"  <xs:element name=\"root\">"
							"    <xs:complexType>"
							"      <xs:sequence>"
							"        <xs:element name=\"child\" type=\"xs:int\"/>"
							"      </xs:sequence>"
							"    </xs:complexType>"
							"  </xs:element>"
							"</xs:schema>";
	const auto xsd_bytes = toBytes(xsd);
	REQUIRE(doc.validateWithXsd(xsd_bytes));
}

TEST_CASE("XmlDocument fails XSD validation for invalid value", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><child>not-an-int</child></root>");
	REQUIRE(doc.loadString(xml));

	const std::string xsd = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
							"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">"
							"  <xs:element name=\"root\">"
							"    <xs:complexType>"
							"      <xs:sequence>"
							"        <xs:element name=\"child\" type=\"xs:int\"/>"
							"      </xs:sequence>"
							"    </xs:complexType>"
							"  </xs:element>"
							"</xs:schema>";
	const auto xsd_bytes = toBytes(xsd);
	REQUIRE_THROWS_AS(doc.validateWithXsd(xsd_bytes), XmlException);
}

TEST_CASE("XmlDocument fails with malformed XSD", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = buildXml("<root><child>9</child></root>");
	REQUIRE(doc.loadString(xml));

	const std::string bad_xsd = "<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">";
	const auto xsd_bytes = toBytes(bad_xsd);
	REQUIRE_THROWS_AS(doc.validateWithXsd(xsd_bytes), XmlException);
}

TEST_CASE("XmlDocument loadString fails on malformed XML", "[serial][xml]")
{
	XmlDocument doc;
	const std::string xml = "<root><child></root>";
	REQUIRE_FALSE(doc.loadString(xml));
	REQUIRE(doc.dumpToString().empty());
}

TEST_CASE("XmlDocument loadFile fails for missing file", "[serial][xml]")
{
	XmlDocument doc;
	REQUIRE_FALSE(doc.loadFile("missing_libxml_wrapper.xml"));
	REQUIRE_FALSE(doc.saveFile("should_not_exist.xml"));
}

TEST_CASE("mergeXmlDocuments merges element children", "[serial][xml]")
{
	XmlDocument main_doc;
	XmlDocument included_doc;

	const XmlElement main_root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	main_root.addChild("main").setText("1");
	main_doc.setRootElement(main_root);

	const XmlElement included_root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	included_root.addChild("a").setText("2");
	included_root.addChild("b").setText("3");
	included_doc.setRootElement(included_root);

	mergeXmlDocuments(main_doc, included_doc);
	const XmlElement merged_root = main_doc.getRootElement();

	REQUIRE(merged_root.childElement("main", 0).getText() == "1");
	REQUIRE(merged_root.childElement("a", 0).getText() == "2");
	REQUIRE(merged_root.childElement("b", 0).getText() == "3");
}

TEST_CASE("removeIncludeElements removes all include nodes", "[serial][xml]")
{
	XmlDocument doc;
	const XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	root.addChild("include").setText("first");
	root.addChild("keep").setText("value");
	root.addChild("include").setText("second");
	doc.setRootElement(root);

	removeIncludeElements(doc);
	const XmlElement reloaded = doc.getRootElement();

	REQUIRE_FALSE(reloaded.childElement("include", 0).isValid());
	REQUIRE(reloaded.childElement("keep", 0).getText() == "value");
}

// TODO: Exercise error paths where libxml2 fails to allocate validation contexts,
// returns a null buffer from xmlDocDumpFormatMemory, or fails to create a new
// xml document in XmlDocument() constructor. These failures are not reliably
// reproducible in unit tests without fault injection.
