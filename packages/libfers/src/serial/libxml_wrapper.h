// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file libxml_wrapper.h
 * @brief Wrapper for managing XML documents and elements using libxml2.
 *
 * This header file provides classes and functions to simplify handling XML documents and elements
 * using the libxml2 library. It includes basic functionality for manipulating XML nodes, attributes,
 * content, and validation using DTD and XSD schemas.
 */

#pragma once

#include <iostream>
#include <libxml/parser.h>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "core/logging.h"
#include "libxml/globals.h"
#include "libxml/xmlstring.h"

/**
 * @class XmlException
 * @brief Exception class for handling XML-related errors.
 */
class XmlException final : public std::runtime_error
{
public:
	/**
	 * @brief Constructor for XmlException.
	 *
	 * @param message The error message associated with the exception.
	 */
	explicit XmlException(const std::string_view message) : std::runtime_error(std::string(message)) {}
};

namespace xml_detail
{
	class XmlCharBuffer
	{
		std::vector<xmlChar> _value;

	public:
		explicit XmlCharBuffer(const std::string_view value)
		{
			_value.reserve(value.size() + 1);
			for (const char character : value)
			{
				_value.push_back(static_cast<xmlChar>(static_cast<unsigned char>(character)));
			}
			_value.push_back(0);
		}

		[[nodiscard]] const xmlChar* c_str() const noexcept { return _value.data(); }

		[[nodiscard]] int contentLength() const
		{
			const auto length = _value.empty() ? static_cast<std::size_t>(0) : _value.size() - 1;
			if (length > static_cast<std::size_t>(std::numeric_limits<int>::max()))
			{
				throw XmlException("XML string exceeds libxml2 length limit.");
			}
			return static_cast<int>(length);
		}
	};

	[[nodiscard]] inline std::string toString(const xmlChar* value)
	{
		if (value == nullptr)
		{
			return "";
		}

		std::string result;
		for (const xmlChar* cursor = value; *cursor != 0; ++cursor)
		{
			result.push_back(static_cast<char>(*cursor));
		}
		return result;
	}

	[[nodiscard]] inline std::string toString(const xmlChar* value, const int length)
	{
		if (value == nullptr || length <= 0)
		{
			return "";
		}

		std::string result;
		result.reserve(static_cast<std::size_t>(length));
		for (int index = 0; index < length; ++index)
		{
			result.push_back(static_cast<char>(value[index]));
		}
		return result;
	}

	[[nodiscard]] inline bool equals(const xmlChar* xml_value, const std::string_view value)
	{
		if (xml_value == nullptr)
		{
			return false;
		}

		for (std::size_t index = 0; index < value.size(); ++index)
		{
			if (xml_value[index] == 0 || static_cast<char>(xml_value[index]) != value[index])
			{
				return false;
			}
		}
		return xml_value[value.size()] == 0;
	}

	[[nodiscard]] inline xmlNodePtr createNode(const std::string_view name)
	{
		const XmlCharBuffer xml_name(name);
		xmlNodePtr node = xmlNewNode(nullptr, xml_name.c_str());
		if (node == nullptr)
		{
			throw XmlException("Failed to create XML node: " + std::string(name));
		}
		return node;
	}

	[[nodiscard]] inline xmlDocPtr createDocument()
	{
		const XmlCharBuffer version("1.0");
		return xmlNewDoc(version.c_str());
	}
}

/**
 * @class XmlElement
 * @brief Class representing a node in an XML document.
 *
 * This class encapsulates an XML element, allowing users to access and manipulate
 * element names, attributes, content, and children. It uses libxml2 for all operations
 * and provides simplified methods to interact with XML nodes.
 */
class XmlElement
{
	xmlNodePtr _node; ///< Pointer to the XML node represented by this object.

public:
	/**
	 * @brief Constructor for XmlElement.
	 *
	 * @param node The xmlNode pointer representing the XML element.
	 */
	explicit XmlElement(xmlNodePtr node) : _node(node) {}

	XmlElement(const XmlElement&) = default;

	XmlElement(XmlElement&&) noexcept = default;

	XmlElement& operator=(const XmlElement&) = default;

	XmlElement& operator=(XmlElement&&) noexcept = default;

	~XmlElement() = default;

	/**
	 * @brief Get the name of the XML element.
	 *
	 * @return The name of the XML element.
	 */
	[[nodiscard]] std::string name() const { return _node == nullptr ? "" : xml_detail::toString(_node->name); }

	/**
	 * @brief Create a new XML element by name.
	 *
	 * @param name The name of the XML element.
	 * @return The newly created XML element.
	 */
	[[nodiscard]] static XmlElement create(const std::string_view name)
	{
		return XmlElement(xml_detail::createNode(name));
	}

	/**
	 * @brief Get the text content of the XML element.
	 *
	 * @return The text content as a string.
	 */
	[[nodiscard]] std::string getText() const
	{
		if (_node == nullptr)
		{
			return "";
		}
		xmlChar* text = xmlNodeGetContent(_node);
		std::string result = xml_detail::toString(text);
		xmlFree(text);
		return result;
	}

	/**
	 * @brief Set the text content of the XML element.
	 *
	 * @param text The text to set as the content of the node.
	 */
	void setText(const std::string_view text) const
	{
		const xml_detail::XmlCharBuffer xml_text(text);
		xmlNodeSetContentLen(_node, xml_text.c_str(), xml_text.contentLength());
	}

	/**
	 * @brief Get the value of an attribute safely.
	 *
	 * @param element The XmlElement to retrieve the attribute from.
	 * @param name The name of the attribute.
	 * @return The value of the attribute.
	 * @throws XmlException if the attribute is not found.
	 */
	static std::string getSafeAttribute(const XmlElement& element, const std::string_view name)
	{
		std::string value;
		const xml_detail::XmlCharBuffer xml_name(name);
		if (xmlChar* attr = xmlGetProp(element.getNode(), xml_name.c_str()))
		{
			value = xml_detail::toString(attr);
			xmlFree(attr);
		}
		else
		{
			throw XmlException("Attribute not found: " + std::string(name));
		}
		return value;
	}

	/**
	 * @brief Get the value of an optional attribute.
	 *
	 * @param element The XmlElement to retrieve the attribute from.
	 * @param name The name of the attribute.
	 * @return The attribute value if present.
	 */
	static std::optional<std::string> getOptionalAttribute(const XmlElement& element, const std::string_view name)
	{
		if (!element.isValid())
		{
			return std::nullopt;
		}
		const xml_detail::XmlCharBuffer xml_name(name);
		if (xmlChar* attr = xmlGetProp(element.getNode(), xml_name.c_str()))
		{
			std::string value = xml_detail::toString(attr);
			xmlFree(attr);
			return value;
		}
		return std::nullopt;
	}

	/**
	 * @brief Set an attribute on the XML element.
	 *
	 * @param name The name of the attribute.
	 * @param value The value to set for the attribute.
	 */
	void setAttribute(const std::string_view name, const std::string_view value) const
	{
		const xml_detail::XmlCharBuffer xml_name(name);
		const xml_detail::XmlCharBuffer xml_value(value);
		xmlSetProp(_node, xml_name.c_str(), xml_value.c_str());
	}

	/**
	 * @brief Add a child element to the current node.
	 *
	 * @param name The name of the new child element.
	 * @return The newly created XmlElement.
	 */
	[[nodiscard]] XmlElement addChild(const std::string_view name) const
	{
		xmlNodePtr child = xml_detail::createNode(name);
		xmlAddChild(_node, child);
		return XmlElement(child);
	}

	/**
	 * @brief Retrieve a child element by name and index.
	 *
	 * @param name The name of the child element (optional).
	 * @param index The index of the child to retrieve.
	 * @return The child element or an invalid XmlElement if not found.
	 */
	[[nodiscard]] XmlElement childElement(const std::string_view name = "", const unsigned index = 0) const noexcept
	{
		if (_node == nullptr)
		{
			return XmlElement(nullptr);
		}
		unsigned count = 0;
		for (auto* child = _node->children; child != nullptr; child = child->next)
		{
			if (child->type == XML_ELEMENT_NODE && (name.empty() || xml_detail::equals(child->name, name)))
			{
				if (count == index)
				{
					return XmlElement(child);
				}
				++count;
			}
		}
		return XmlElement(nullptr);
	}

	/**
	 * @brief Check if the XML element is valid.
	 *
	 * @return True if the node is valid, otherwise false.
	 */
	[[nodiscard]] bool isValid() const noexcept { return _node != nullptr; }

	/**
	 * @brief Get the underlying XML node pointer.
	 *
	 * @return The underlying XML node pointer.
	 */
	[[nodiscard]] xmlNodePtr getNode() const noexcept { return _node; }
};

/**
 * @class XmlDocument
 * @brief Class for managing XML documents.
 */
class XmlDocument
{
	std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> _doc; ///< Pointer to the XML document.

public:
	/**
	 * @brief Constructor for XmlDocument.
	 *
	 * @throws std::runtime_error if the document creation fails.
	 */
	XmlDocument() : _doc(xml_detail::createDocument(), &xmlFreeDoc)
	{
		if (!_doc)
		{
			throw XmlException("Failed to create XML document.");
		}
	}

	~XmlDocument() = default;

	XmlDocument(const XmlDocument&) = delete;

	XmlDocument(XmlDocument&&) noexcept = default;

	XmlDocument& operator=(const XmlDocument&) = delete;

	XmlDocument& operator=(XmlDocument&&) noexcept = default;

	/**
	 * @brief Load an XML file into the document.
	 *
	 * @param filename The name of the file to load.
	 * @return True if the file was successfully loaded, otherwise false.
	 */
	[[nodiscard]] bool loadFile(std::string_view filename);

	/**
	 * @brief Load an XML document from a string in memory.
	 *
	 * @param content The string containing the XML document.
	 * @return True if the string was successfully parsed, otherwise false.
	 */
	[[nodiscard]] bool loadString(const std::string& content);

	/**
	 * @brief Save the document to a file.
	 *
	 * @param filename The name of the file to save to.
	 * @return True if the file was successfully saved, otherwise false.
	 */
	[[nodiscard]] bool saveFile(const std::string_view filename) const
	{
		if (!_doc)
		{
			LOG(logging::Level::ERROR, "Document is null; Cannot save to file");
			return false;
		}
		return xmlSaveFormatFileEnc(filename.data(), _doc.get(), "UTF-8", 1) != -1;
	}

	/**
	 * @brief Dumps the document to a string.
	 *
	 * @return A string containing the XML document.
	 */
	[[nodiscard]] std::string dumpToString() const;

	/**
	 * @brief Set the root element of the document.
	 *
	 * @param root The root element to set.
	 * @throws std::runtime_error if the document is not created.
	 */
	void setRootElement(const XmlElement& root) const
	{
		if (!_doc)
		{
			throw std::runtime_error("Document not created");
		}
		xmlDocSetRootElement(_doc.get(), root.getNode());
	}

	/**
	 * @brief Get the root element of the document.
	 *
	 * @return The root element.
	 * @throws std::runtime_error if the document is not loaded or the root element is missing.
	 */
	[[nodiscard]] XmlElement getRootElement() const
	{
		if (!_doc)
		{
			throw std::runtime_error("Document not loaded");
		}
		xmlNodePtr root = xmlDocGetRootElement(_doc.get());
		if (root == nullptr)
		{
			throw std::runtime_error("Root element not found");
		}
		return XmlElement(root);
	}

	/**
	 * @brief Validate the document using a DTD.
	 *
	 * @param dtdData The DTD data used for validation.
	 * @return True if the document is valid according to the DTD.
	 * @throws XmlException if the DTD is invalid or the validation fails.
	 */
	[[nodiscard]] bool validateWithDtd(std::span<const unsigned char> dtdData) const;

	/**
	 * @brief Validate the document using an XSD schema.
	 *
	 * @param xsdData The XSD data used for validation.
	 * @return True if the document is valid according to the XSD schema.
	 * @throws XmlException if the XSD is invalid or the validation fails.
	 */
	[[nodiscard]] bool validateWithXsd(std::span<const unsigned char> xsdData) const;
};

/**
 * @brief Merge two XML documents.
 *
 * @param mainDoc The main XML document.
 * @param includedDoc The XML document to include.
 */
void mergeXmlDocuments(const XmlDocument& mainDoc, const XmlDocument& includedDoc);

/**
 * @brief Remove "include" elements from the XML document.
 *
 * @param doc The XML document from which to remove the "include" elements.
 */
void removeIncludeElements(const XmlDocument& doc);
