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

#include <cctype>
#include <cstdarg>
#include <format>
#include <string>

#include "libxml/encoding.h"
#include "libxml/parser.h"
#include "libxml/valid.h"
#include "libxml/xmlIO.h"
#include "libxml/xmlerror.h"
#include "libxml/xmlschemas.h"

namespace
{
	/// Captures libxml2 generic errors into a std::string buffer.
	void libxmlGenericErrorCallback(void* ctx, const char* msg, ...)
	{
		auto* err_str = static_cast<std::string*>(ctx);
		if ((err_str == nullptr) || (msg == nullptr))
		{
			return;
		}

		char buf[1024];
		va_list args;
		va_start(args, msg);
		vsnprintf(buf, sizeof(buf), msg, args);
		va_end(args);

		err_str->append(buf);
	}

	/// Formats a raw XML error message with a title and remediation hint.
	std::string formatError(const std::string_view title, const std::string& rawErrors, const std::string_view hint)
	{
		std::string errors = rawErrors;

		// Strip trailing newlines or whitespace
		while (!errors.empty() && (std::isspace(static_cast<unsigned char>(errors.back())) != 0))
		{
			errors.pop_back();
		}

		// Indent multiline errors to align nicely within the box
		size_t pos = 0;
		while ((pos = errors.find('\n', pos)) != std::string::npos)
		{
			errors.replace(pos, 1, "\n│      ");
			pos += 8;
		}

		return std::format("\n"
						   "┌─ {} \n"
						   "│\n"
						   "│ Error Details:\n"
						   "│      {}\n"
						   "│\n"
						   "│ Hint: {}\n"
						   "└────────────────────────────────────────────────────────────────────────────",
						   title, errors.empty() ? "Unknown validation error." : errors, hint);
	}

	/// Retrieves the most recent libxml parser error as a formatted string.
	std::string getXmlLastErrorFormatted()
	{
		const xmlError* err = xmlGetLastError();
		if ((err != nullptr) && (err->message != nullptr))
		{
			std::string msg = err->message;
			while (!msg.empty() && (std::isspace(static_cast<unsigned char>(msg.back())) != 0))
			{
				msg.pop_back();
			}
			if (err->line > 0)
			{
				return std::format("Line {}: {}", err->line, msg);
			}
			return msg;
		}
		return "Syntax error or malformed XML.";
	}
}

bool XmlDocument::validateWithDtd(const std::span<const unsigned char> dtdData) const
{
	xmlDtdPtr dtd =
		xmlIOParseDTD(nullptr,
					  xmlParserInputBufferCreateMem(reinterpret_cast<const char*>(dtdData.data()),
													static_cast<int>(dtdData.size()), XML_CHAR_ENCODING_UTF8),
					  XML_CHAR_ENCODING_UTF8);
	if (dtd == nullptr)
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

	// Bind our custom error handler into the DTD Validation Context
	std::string dtdErrors;
	validation_ctxt->userData = &dtdErrors;
	validation_ctxt->error = libxmlGenericErrorCallback;
	validation_ctxt->warning = libxmlGenericErrorCallback;

	const bool is_valid = xmlValidateDtd(validation_ctxt.get(), _doc.get(), dtd) != 0;
	xmlFreeDtd(dtd);

	if (!is_valid)
	{
		std::string fancyError = formatError("XML DTD Validation Failed", dtdErrors,
											 "Check your scenario XML tags and attributes against 'fers-xml.dtd'.");
		LOG(logging::Level::ERROR, "{}", fancyError);
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

	// Bind custom error handler into the Schema Parse Context
	std::string xsdParseErrors;
	xmlSchemaSetParserErrors(schema_parser_ctxt.get(), libxmlGenericErrorCallback, libxmlGenericErrorCallback,
							 &xsdParseErrors);

	const std::unique_ptr<xmlSchema, decltype(&xmlSchemaFree)> schema(xmlSchemaParse(schema_parser_ctxt.get()),
																	  xmlSchemaFree);
	if (!schema)
	{
		std::string fancyError =
			formatError("XSD Schema Parse Failed", xsdParseErrors, "The internal XSD schema is invalid.");
		LOG(logging::Level::FATAL, "{}", fancyError);
		throw XmlException("Failed to parse schema from memory.");
	}

	const std::unique_ptr<xmlSchemaValidCtxt, decltype(&xmlSchemaFreeValidCtxt)> schema_valid_ctxt(
		xmlSchemaNewValidCtxt(schema.get()), xmlSchemaFreeValidCtxt);
	if (!schema_valid_ctxt)
	{
		throw XmlException("Failed to create schema validation context.");
	}

	// Bind custom error handler into the Schema Validation Context
	std::string xsdErrors;
	xmlSchemaSetValidErrors(schema_valid_ctxt.get(), libxmlGenericErrorCallback, libxmlGenericErrorCallback,
							&xsdErrors);

	if (const bool is_valid = xmlSchemaValidateDoc(schema_valid_ctxt.get(), _doc.get()) == 0; !is_valid)
	{
		std::string fancyError = formatError("XML XSD Validation Failed", xsdErrors,
											 "Check your scenario XML tags and attributes against 'fers-xml.xsd'.");
		LOG(logging::Level::ERROR, "{}", fancyError);
		throw XmlException("XML failed XSD validation.");
	}

	return true;
}

void mergeXmlDocuments(const XmlDocument& mainDoc, const XmlDocument& includedDoc)
{
	const XmlElement main_root = mainDoc.getRootElement();
	const XmlElement included_root = includedDoc.getRootElement();

	for (xmlNodePtr child = included_root.getNode()->children; child != nullptr; child = child->next)
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
	xmlResetLastError();
	// Pass NOERROR and NOWARNING to prevent default terminal spam, so we handle it cleanly
	_doc.reset(xmlReadFile(filename.data(), nullptr, XML_PARSE_NOERROR | XML_PARSE_NOWARNING));
	if (!_doc)
	{
		std::string fancyError =
			formatError("XML Parsing Failed", getXmlLastErrorFormatted(), "Ensure the XML file is well-formed.");
		LOG(logging::Level::ERROR, "{}", fancyError);
		return false;
	}
	return true;
}

bool XmlDocument::loadString(const std::string& content)
{
	xmlResetLastError();
	_doc.reset(xmlReadMemory(content.c_str(), static_cast<int>(content.length()), "in_memory.xml", nullptr,
							 XML_PARSE_NOERROR | XML_PARSE_NOWARNING));
	if (!_doc)
	{
		std::string fancyError =
			formatError("XML Parsing Failed", getXmlLastErrorFormatted(), "Ensure the XML string is well-formed.");
		LOG(logging::Level::ERROR, "{}", fancyError);
		return false;
	}
	return true;
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
	if (buffer == nullptr)
	{
		LOG(logging::Level::ERROR, "Failed to dump XML document to memory buffer.");
		return "";
	}
	const std::string result(reinterpret_cast<const char*>(buffer), static_cast<size_t>(size));
	xmlFree(buffer);
	return result;
}
