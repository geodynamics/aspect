# Copyright (C) 2022 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASPECT; see the file doc/COPYING.  If not see
# <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

MESSAGE(STATUS "===== Configuring ASPECT documentation =============")

FIND_PACKAGE(Perl)


####################################################
# Generate input files for the manual
####################################################

# Step 1: Find all of the .prm files that will be used in
#         the manual. Do the same for the .tex and .md
#         files
FILE(GLOB _doc_directories_cookbooks
     ${CMAKE_SOURCE_DIR}/cookbooks/*/doc)
FILE(GLOB _doc_directories_benchmarks
     ${CMAKE_SOURCE_DIR}/benchmarks/*/doc)
SET(_doc_directories
    ${CMAKE_SOURCE_DIR}/doc/manual
    ${CMAKE_SOURCE_DIR}/doc/logo)

LIST(APPEND _doc_directories
    ${_doc_directories_cookbooks}
    ${_doc_directories_benchmarks})

# Collect the .tex, .md, and .prm files
FOREACH(_dir ${_doc_directories})
  FILE(GLOB_RECURSE _my_dot_tex_files
       RELATIVE ${CMAKE_SOURCE_DIR}
       ${_dir}/*tex)
  LIST(APPEND _dot_tex_files
       ${_my_dot_tex_files})

  FILE(GLOB_RECURSE _my_dot_md_files
       RELATIVE ${CMAKE_SOURCE_DIR}
       ${_dir}/*md)
  LIST(APPEND _dot_md_files
       ${_my_dot_md_files})

  FILE(GLOB_RECURSE _my_dot_prm_files
       RELATIVE ${CMAKE_SOURCE_DIR}
       ${_dir}/*prm)
  LIST(APPEND _dot_prm_files
       ${_my_dot_prm_files})
ENDFOREACH()


# Step 2: Now also find all of the image files that will find their
#         way into the manual. These can have file endings .png
#         or .svg.
# Collect the .png and .svg files
FOREACH(_dir ${_doc_directories})
  FILE(GLOB_RECURSE _my_image_files
       RELATIVE ${CMAKE_SOURCE_DIR}
       ${_dir}/*png
       ${_dir}/*svg)
  LIST(APPEND _image_files
       ${_my_image_files})
ENDFOREACH()


# Step 3: Set up targets to run each of the .prm files
#         through perl and process them
ADD_CUSTOM_TARGET(build_dot_prm_dot_tex)

FOREACH(_dot_prm ${_dot_prm_files})
  # Set up a command to generate the .prm.tex file
  #
  # We do this using the annotate.pl script, and then processing
  # index entries to contain at most three levels using a
  # second perl script.
  ADD_CUSTOM_COMMAND(
    OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    COMMAND
      ${PERL_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/annotate.pl
          ${CMAKE_SOURCE_DIR}/${_dot_prm}
          > ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    COMMAND
      ${PERL_EXECUTABLE}
         -pi
          ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/escape_index_entries.pl
          ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    DEPENDS
      ${CMAKE_SOURCE_DIR}/${_dot_prm}
      ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/annotate.pl
      ${CMAKE_CURRENT_SOURCE_DIR}/manual/cookbooks/escape_index_entries.pl
    COMMENT
      "Building marked up version of ${_dot_prm} in latex format."
    )

  # Then we also need a target for the operation above. Target
  # names may not contain slashes, so we need to escape those
  STRING(REGEX REPLACE "[^a-zA-Z0-9]" "" _dot_prm_escaped "${_dot_prm}")
  ADD_CUSTOM_TARGET(build_dot_prm_dot_tex_${_dot_prm_escaped}
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${_dot_prm}.tex
    )

  # Finally add the target to the dependencies of the higher-level
  # target
  ADD_DEPENDENCIES(build_dot_prm_dot_tex build_dot_prm_dot_tex_${_dot_prm_escaped})
ENDFOREACH()


# Step 4: Get all of the latex source files into the current directory. There
#         may be a way to get around doing this and just sourcing the latex
#         sources from the source directory while getting the .prm.tex files
#         from the build directory, but that seems unnecessarily complicated.
#         Furthermore, it allows us to do a minimal amount of processing
#         when copying files from the source to the build directory.
ADD_CUSTOM_TARGET(build_dot_tex)

FOREACH(_dot_tex ${_dot_tex_files})
  # Set up a command to generate the .tex file from the source .tex file.
  #
  # Right now, this essentially involves copying the source file to the build
  # directory. The only modification we make are:
  # * historically, we have called the processed .prm files .prm.out, whereas
  #   above we build them as .prm.tex. Make the appropriate modification in
  #   the .tex files.
  # * When we include these .prm.tex files from the cookbooks, we do this
  #   through symlinks the git repo contains in doc/manual/cookbooks, but
  #   these are not copied into the build directory. So we need to provide
  #   correct paths to these sorts of places. The exception is the path
  #   doc/manual/cookbooks/overview/, which is a real directory.
  ADD_CUSTOM_COMMAND(
    OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/${_dot_tex}
    COMMAND
      cat ${CMAKE_SOURCE_DIR}/${_dot_tex}
        | ${PERL_EXECUTABLE} -p -e 's/\.prm\.out\}/.prm.tex\}/g;'
        > ${CMAKE_CURRENT_BINARY_DIR}/${_dot_tex}
    DEPENDS
      ${CMAKE_SOURCE_DIR}/${_dot_tex}
    COMMENT
      "Copying ${_dot_tex} to build directory."
    )

  # Then we also need a target for the operation above. Target
  # names may not contain slashes, so we need to escape those
  STRING(REGEX REPLACE "[^a-zA-Z0-9]" "" _dot_tex_escaped "${_dot_tex}")
  ADD_CUSTOM_TARGET(build_dot_tex_${_dot_tex_escaped}
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${_dot_tex}
    )

  # Finally add the target to the dependencies of the higher-level
  # target
  ADD_DEPENDENCIES(build_dot_tex build_dot_tex_${_dot_tex_escaped})
ENDFOREACH()


# Step 5: Do the same with image files
ADD_CUSTOM_TARGET(build_images)

FOREACH(_image ${_image_files})
  # Set up a command to generate the .tex file from the source .tex file.
  #
  # Right now, this simply involves copying the source file to the build
  # directory.
  ADD_CUSTOM_COMMAND(
    OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/${_image}
    COMMAND
      cp
          ${CMAKE_SOURCE_DIR}/${_image}
          ${CMAKE_CURRENT_BINARY_DIR}/${_image}
    DEPENDS
      ${CMAKE_SOURCE_DIR}/${_image}
    COMMENT
      "Copying ${_image} to build directory."
    )

  # Then we also need a target for the operation above. Target
  # names may not contain slashes, so we need to escape those
  STRING(REGEX REPLACE "[^a-zA-Z0-9]" "" _image_escaped "${_image}")
  ADD_CUSTOM_TARGET(build_image_${_image_escaped}
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${_image}
    )

  # Finally add the target to the dependencies of the higher-level
  # target
  ADD_DEPENDENCIES(build_images build_image_${_image_escaped})
ENDFOREACH()


# Step 6: We have one image file that we need to create -- the plugin
#         graph. This happens by calling the 'neato' program with
#         the plugin graph input file generated by ASPECT.
#
#         Right now, this file is put into the repository, so we don't
#         need to generate it. As a consequence, just add the rule
#         whenever actually necessary. (If the file is found, then it's
#         already copied into the build directory in step 5 above.)
IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/manual/plugin_graph.svg)
  ADD_CUSTOM_COMMAND(
    OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/plugin_graph.svg
    COMMAND
      neato
        ${CMAKE_CURRENT_SOURCE_DIR}/manual/plugin_graph.dot
        -Tsvg
        -o ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/plugin_graph.svg
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/manual/plugin_graph.dot
    COMMENT
      "Generating the plugin graph image."
  )

  ADD_CUSTOM_TARGET(build_plugin_graph_svg
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/plugin_graph.svg
  )

  # Add the step above as a dependency -- somewhat randomly -- to the
  # image generation target above. This way we don't have to do another
  # IF(NOT EXISTS...) when adding dependencies to the larger targets below.
  ADD_DEPENDENCIES(build_images
                   build_plugin_graph_svg)
ENDIF()


# Step 7: Another preparatory step: Copy the manual.bib file into the
#         directory where we will build the manual.
ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/manual.bib
  COMMAND
    cp
      ${CMAKE_CURRENT_SOURCE_DIR}/manual/manual.bib
      ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/manual.bib
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/manual/manual.bib
    COMMENT
      "Copying the bibliography database."
  )



# Step 8: Run pdflatex to generate the pdf version of the manual. This
#         obviously depends on all of the previous targets.
ADD_CUSTOM_COMMAND(
    OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/manual.pdf
    WORKING_DIRECTORY
      doc/manual
    COMMAND
      rm -r svg-inkscape
    COMMAND
      pdflatex --shell-escape --interaction=batchmode manual.tex
    COMMAND
      bibtex manual
    COMMAND
      makeindex prmindex
    COMMAND
      pdflatex --shell-escape --interaction=batchmode manual.tex
    COMMAND
      pdflatex --shell-escape --interaction=batchmode manual.tex
    COMMAND
      pdflatex --shell-escape --interaction=batchmode manual.tex
    DEPENDS
      build_dot_tex
      build_images
      build_dot_prm_dot_tex
      ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/manual.bib
    COMMENT
      "Generating manual.pdf."
  )

ADD_CUSTOM_TARGET(build_manual_tex
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/doc/manual/manual.pdf
)





################################################################
#
# The final target: Build it all
#
################################################################
ADD_CUSTOM_TARGET(documentation
  DEPENDS
    build_manual_tex
)
