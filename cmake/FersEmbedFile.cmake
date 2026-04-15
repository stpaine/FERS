if (NOT DEFINED INPUT)
	message(FATAL_ERROR "FersEmbedFile requires INPUT")
endif ()

if (NOT DEFINED OUTPUT)
	message(FATAL_ERROR "FersEmbedFile requires OUTPUT")
endif ()

if (NOT DEFINED VAR)
	message(FATAL_ERROR "FersEmbedFile requires VAR")
endif ()

file(READ "${INPUT}" _fers_embed_hex HEX)
string(LENGTH "${_fers_embed_hex}" _fers_embed_hex_length)
math(EXPR _fers_embed_byte_count "${_fers_embed_hex_length} / 2")

set(_fers_embed_body)
if (_fers_embed_byte_count GREATER 0)
	math(EXPR _fers_embed_last_index "${_fers_embed_byte_count} - 1")
	foreach (_fers_embed_index RANGE 0 ${_fers_embed_last_index})
		math(EXPR _fers_embed_pos "${_fers_embed_index} * 2")
		string(SUBSTRING "${_fers_embed_hex}" ${_fers_embed_pos} 2 _fers_embed_byte)

		if (_fers_embed_index GREATER 0)
			string(APPEND _fers_embed_body ",")
		endif ()

		math(EXPR _fers_embed_col "${_fers_embed_index} % 12")
		if (_fers_embed_col EQUAL 0)
			string(APPEND _fers_embed_body "\n\t")
		else ()
			string(APPEND _fers_embed_body " ")
		endif ()

		string(APPEND _fers_embed_body "0x${_fers_embed_byte}")
	endforeach ()
	string(APPEND _fers_embed_body "\n")
endif ()

file(WRITE "${OUTPUT}"
	 "unsigned char ${VAR}[] = {${_fers_embed_body}};\n"
	 "unsigned int ${VAR}_len = ${_fers_embed_byte_count};\n"
)
