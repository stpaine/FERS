// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file libxml_wrapper.cpp
 * @brief Wrapper for managing XML documents and elements using libxml2.
 */

#include "libxml_wrapper.h"

#include "libxml/encoding.h"
#include "libxml/parser.h"
#include "libxml/valid.h"
#include "libxml/xmlIO.h"
#include "libxml/xmlschemas.h"

bool XmlDocument::validateWithDtd(const std::span<const unsigned char> dtdData) const
{
	xmlDtdPtr dtd =
		xmlIOParseDTD(nullptr,
					  xmlParserInputBufferCreateMem(reinterpret_cast<const char*>(dtdData.data()),
													static_cast<int>(dtdData.size()), XML_CHAR_ENCODING_UTF8),
					  XML_CHAR_ENCODING_UTF8);
	if (!dtd)
	{
		throw XmlException("Failed to parse DTD from memory.");
	}

	const std::unique_ptr<xmlValidCtxt, decltype(&xmlFreeValidCtxt)> validation_ctxt(xmlNewValidCtxt(),
																					 xmlFreeValidCtxt);
	if (!validation_ctxt)
	{
		xmlFreeDtd(dtd);
		throw XmlException("Failed to create validation context.");
	}

	const bool is_valid = xmlValidateDtd(validation_ctxt.get(), _doc.get(), dtd);
	xmlFreeDtd(dtd);

	if (!is_valid)
	{
		throw XmlException("XML failed DTD validation.");
	}

	return true;
}

bool XmlDocument::validateWithXsd(const std::span<const unsigned char> xsdData) const
{
	const std::unique_ptr<xmlSchemaParserCtxt, decltype(&xmlSchemaFreeParserCtxt)> schema_parser_ctxt(
		xmlSchemaNewMemParserCtxt(reinterpret_cast<const char*>(xsdData.data()), static_cast<int>(xsdData.size())),
		xmlSchemaFreeParserCtxt);
	if (!schema_parser_ctxt)
	{
		throw XmlException("Failed to create schema parser context.");
	}

	const std::unique_ptr<xmlSchema, decltype(&xmlSchemaFree)> schema(xmlSchemaParse(schema_parser_ctxt.get()),
																	  xmlSchemaFree);
	if (!schema)
	{
		throw XmlException("Failed to parse schema from memory.");
	}

	const std::unique_ptr<xmlSchemaValidCtxt, decltype(&xmlSchemaFreeValidCtxt)> schema_valid_ctxt(
		xmlSchemaNewValidCtxt(schema.get()), xmlSchemaFreeValidCtxt);
	if (!schema_valid_ctxt)
	{
		throw XmlException("Failed to create schema validation context.");
	}

	if (const bool is_valid = xmlSchemaValidateDoc(schema_valid_ctxt.get(), _doc.get()) == 0; !is_valid)
	{
		throw XmlException("XML failed XSD validation.");
	}

	return true;
}

void mergeXmlDocuments(const XmlDocument& mainDoc, const XmlDocument& includedDoc)
{
	const XmlElement main_root = mainDoc.getRootElement();
	const XmlElement included_root = includedDoc.getRootElement();

	for (xmlNodePtr child = included_root.getNode()->children; child; child = child->next)
	{
		if (child->type == XML_ELEMENT_NODE)
		{
			xmlNodePtr new_node = xmlCopyNode(child, 1);
			xmlAddChild(main_root.getNode(), new_node);
		}
	}
}

void removeIncludeElements(const XmlDocument& doc)
{
	const XmlElement root = doc.getRootElement();

	while (true)
	{
		if (XmlElement include_element = root.childElement("include", 0); include_element.isValid())
		{
			xmlUnlinkNode(include_element.getNode());
			xmlFreeNode(include_element.getNode());
		}
		else
		{
			break;
		}
	}
}

bool XmlDocument::loadFile(const std::string_view filename)
{
	_doc.reset(xmlReadFile(filename.data(), nullptr, 0));
	return _doc != nullptr;
}

bool XmlDocument::loadString(const std::string& content)
{
	_doc.reset(xmlReadMemory(content.c_str(), static_cast<int>(content.length()), "in_memory.xml", nullptr, 0));
	return _doc != nullptr;
}

std::string XmlDocument::dumpToString() const
{
	if (!_doc)
	{
		LOG(logging::Level::ERROR, "Document is null; Cannot dump to string");
		return "";
	}
	xmlChar* buffer = nullptr;
	int size = 0;
	xmlDocDumpFormatMemory(_doc.get(), &buffer, &size, 1);
	if (!buffer)
	{
		LOG(logging::Level::ERROR, "Failed to dump XML document to memory buffer.");
		return "";
	}
	const std::string result(reinterpret_cast<const char*>(buffer), static_cast<size_t>(size));
	xmlFree(buffer);
	return result;
}
