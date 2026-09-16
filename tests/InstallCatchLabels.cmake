foreach(target_name IN LISTS ALL_TEST_TARGETS)
  file(GLOB discovered_file "${BUILD_DIR}/${target_name}-*_tests.cmake")
  if(NOT discovered_file)
    continue()
  endif()

  # catch_discover_tests writes LABELS as a bare word (LABELS foo) when there is
  # a single tag with no characters needing escaping, and otherwise as a CMake
  # bracket argument (LABELS [==[tag1;tag2]==]).
  file(READ "${discovered_file}" discovered_content)
  string(REGEX MATCHALL "LABELS ([-./:a-zA-Z0-9_]+|\\[==\\[[^]]*\\]==\\])" label_props "${discovered_content}")
  set(target_labels "")
  foreach(label_prop ${label_props})
    string(REGEX REPLACE "^LABELS " "" label_list "${label_prop}")
    string(REGEX REPLACE "^\\[==\\[|\\]==\\]$" "" label_list "${label_list}")
    list(APPEND target_labels ${label_list})
  endforeach()

  if(target_labels)
    list(REMOVE_DUPLICATES target_labels)
    file(APPEND "${INSTALLED_CTEST_FILE}"
         "set_tests_properties(${target_name} PROPERTIES LABELS \"${target_labels}\")\n")
  endif()
endforeach()
