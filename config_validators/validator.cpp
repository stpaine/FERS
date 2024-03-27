// FERS input Validator Function
// Makes calls to xml_validator.cpp and kml_visualiser.cpp
// Script written by Michael Altshuler 
// University of Cape Town: ALTMIC003

// The following implementation makes use of the Document Object Model (DOM) and a custom XSD Schema.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// Including necessary Xerces-C++ header files:
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/validators/schema/SchemaValidator.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>

using namespace xercesc;
using namespace std;

// Error handler that catches validation errors and further provides information about the errors.
class MyErrorHandler : public ErrorHandler {
public:
    void warning(const SAXParseException& e) {
        cerr << "Warning: " << XMLString::transcode(e.getMessage()) << " at line " << e.getLineNumber() << ", column " << e.getColumnNumber() << endl;
    }

    void error(const SAXParseException& e) {
        cerr << "Error: " << XMLString::transcode(e.getMessage()) << " at line " << e.getLineNumber() << ", column " << e.getColumnNumber() << endl;
    }

    void fatalError(const SAXParseException& e) {
        cerr << "Fatal Error: " << XMLString::transcode(e.getMessage()) << " at line " << e.getLineNumber() << ", column " << e.getColumnNumber() << endl;
    }

    void resetErrors() {}
};

// Function to extract the file name from a given path
std::string getFileNameFromPath(const std::string &filePath) {
    size_t lastSlash = filePath.find_last_of("/\\");
    size_t lastDot = filePath.find_last_of(".");
    if (lastSlash == std::string::npos) {
        lastSlash = 0;
    } else {
        lastSlash++;
    }
    return filePath.substr(lastSlash, lastDot - lastSlash);
}


int main(int argc, char* argv[]) {

    try {

        // Check if the user provided the mode argument
        std::string mode = "non-verbose"; // Default mode
        std::string file_path;

        if (argc == 2) {
            file_path = argv[1];
        } else if (argc == 3) {
            std::string mode_arg = argv[1];
            if (mode_arg == "-v" || mode_arg == "--verbose") {
                mode = "verbose";
            } else if (mode_arg == "-q" || mode_arg == "--quiet" || mode_arg == "--non-verbose") {
                mode = "non-verbose";
            } else {
                std::cerr << "Invalid mode argument. Please provide '-v' or '--verbose' for verbose mode, or '-q', '--quiet', or '--non-verbose' for non-verbose mode." << std::endl;
                return 1;
            }
            file_path = argv[2];
        } else {
            std::cerr << "Invalid number of arguments. Please provide a file path and an optional mode (-v for verbose)." << std::endl;
            return 1;
        }



        // Initializing Xerces-C++ library:
        xercesc::XMLPlatformUtils::Initialize();

        // A XercesDOMParser object is set along with its features:
        xercesc::XercesDOMParser parser;
        // Enables validation during parsing.
        parser.setValidationScheme(xercesc::XercesDOMParser::Val_Always);
        // Enables validation against an XSD Schema file.
        parser.setDoNamespaces(true);
        // Enables full schema constraint checking during validation.
        parser.setDoSchema(true);

        // Set the custom XSD file for the parser by specifying the XSD file path and XSD file name
        parser.setExternalNoNamespaceSchemaLocation("/Users/michaelaltshuler/Documents/5th Year/EEE4022F:Thesis/FERS Features/FERS Validator/FERS-schema/fers-xml.xsd");

        // error handler instance created and set on the parser.
        MyErrorHandler errorHandler;
        parser.setErrorHandler(&errorHandler);

        // Establish a DOMDocument object and parse the input FERSXML file:
        parser.parse(file_path.c_str());

        cout<< endl;

        if(parser.getErrorCount() == 0) {
                std::cout << "User is in " << mode << " mode." << std::endl;
                cout << endl;

                std::cout << "XML document is valid" << std::endl;

                // Command to run xml_validator_output with mode and file_path as argument
                std::string command = "./xml_validator_output " + file_path + " " + mode;
                // Run the command using system() function
                int result = std::system(command.c_str());

                if (result == -1) {
                    std::cerr << "Failed to run xml_validator_output." << std::endl;
                    return 1;
                }

                std::cout << endl;

                // Prompting user whether they want to run the kml_visualiser
                char run_kml_visualiser = 'n';
                char custom_coordinates;
                double referenceLatitude = -33.9545;
                double referenceLongitude = 18.4563;
                double referenceAltitude = 0;
                std::string outputName;
                std::string input;

                std::cout << "Do you want to run the kml_visualiser program? (y/N): ";
                std::getline(std::cin, input);

                if (!input.empty()) {
                    run_kml_visualiser = input[0];
                }

                if (run_kml_visualiser == 'y' || run_kml_visualiser == 'Y') {
                    std::cout << "Name a KML output file (defualt file name = file name of fersxml): ";
                    std::getline(std::cin, outputName);

                    if (outputName.empty()) {
                        outputName = getFileNameFromPath(file_path);
                    }

                    std::cout << "Do you want to enter custom referenceLatitude, referenceLongitude, and referenceAltitude? (y/N): ";
                    input.clear();
                    std::getline(std::cin, input);
                    custom_coordinates = input.empty() ? 'n' : input[0];
                    std::cout << "Default coordinates set to University of Cape Town: <133.9577° S, 18.4612° E, 0>" << std::endl;

                    if (custom_coordinates == 'y' || custom_coordinates == 'Y') {
                        std::string input;

                        std::cout << "Enter referenceLatitude (default: -33.9545): ";
                        std::cin >> input;
                        std::istringstream iss_latitude(input);
                        if (!(iss_latitude >> referenceLatitude)) {
                            std::cerr << "Invalid input. Please enter a valid number for referenceLatitude." << std::endl;
                            return 1;
                        }

                        std::cout << "Enter referenceLongitude (default: 18.4563): ";
                        std::cin >> input;
                        std::istringstream iss_longitude(input);
                        if (!(iss_longitude >> referenceLongitude)) {
                            std::cerr << "Invalid input. Please enter a valid number for referenceLongitude." << std::endl;
                            return 1;
                        }

                        std::cout << "Enter referenceAltitude (default: 0): ";
                        std::cin >> input;
                        std::istringstream iss_altitude(input);
                        if (!(iss_altitude >> referenceAltitude)) {
                            std::cerr << "Invalid input. Please enter a valid number for referenceAltitude." << std::endl;
                            return 1;
                        }
                    }

                    std::string kml_visualiser_command = "./kml_visualiser";

                    if (custom_coordinates == 'y' || custom_coordinates == 'Y') {
                        kml_visualiser_command += " " + file_path + " " + outputName + ".kml" + " " + std::to_string(referenceLatitude) + " " + std::to_string(referenceLongitude) + " " + std::to_string(referenceAltitude);
                        
                        int kml_visualiser_result = std::system(kml_visualiser_command.c_str());

                        if (kml_visualiser_result == -1) {
                            std::cerr << "Failed to run kml_visualiser." << std::endl;
                            return 1;
                        }
                    } else {
                        kml_visualiser_command += " " + file_path + " " + outputName + ".kml";
                        int kml_visualiser_result = std::system(kml_visualiser_command.c_str());

                        if (kml_visualiser_result == -1) {
                            std::cerr << "Failed to run kml_visualiser." << std::endl;
                            return 1;
                        }
                    }

                    std::cout << outputName + ".kml" + " outputted to current working directory." << std::endl; 

                } else if (run_kml_visualiser != 'n' && run_kml_visualiser != 'N') {
                    std::cerr << "Invalid input. Please enter 'y' for yes or 'n' for no." << std::endl;
                    return 1;
                }

        }

        else {
            std::cout << "FERSXML invalid: Take note of given errors and warnings above." << std::endl;
        }

            // Function call deallocates memory associated with he parsed FERSXML document.
            //parser.resetDocumentPool();
            // Shuts down and deallocates any memory associated with the Xerces-C++ library.
            //xercesc::XMLPlatformUtils::Terminate();
            return 0;

    } catch (const xercesc::XMLException& e) {
        std::cerr << "Error parsing XML document: " << e.getMessage() << std::endl;
        return 1;

    }  catch (const xercesc::DOMException& e) {
        std::cerr << "Error parsing XML document: " << e.getMessage() << std::endl;
        return 1;

    } catch (const xercesc::SAXParseException& e) {
        std::cerr << "Error parsing XML document: " << e.getMessage() << " at line " << e.getLineNumber() << std::endl;
        return 1;

    } catch (...) {
        std::cerr << "Error parsing XML document: unknown exception" << std::endl;
        return 1;
    }
}