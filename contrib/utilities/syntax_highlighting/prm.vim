" Vim syntax file
" Language:         Deal.II/Aspect input parameter file
" Maintainer:       Jonathan Perry-Houts <jperryh2@uoregon.edu>
" Latest Revision:  2016 03 November

if exists("b:current_syntax")
  finish
endif

"---------------------------------------------------------------------------/
" Comments
"---------------------------------------------------------------------------/
syn keyword prmTodo     containedin=prmComment contained TODO FIXME NOTE
syn match   prmComment  "#.*$" contains=prmTodo

"---------------------------------------------------------------------------/
" Constants
"---------------------------------------------------------------------------/
syn match   prmNumber   "\<0[oO]\=\o\+[Ll]\=\>"
syn match   prmNumber   "\<0[xX]\x\+[Ll]\=\>"
syn match   prmNumber   "\<0[bB][01]\+[Ll]\=\>"
syn match   prmNumber   "\<\%([1-9]\d*\|0\)[Ll]\=\>"
syn match   prmNumber   "\<\d\+[jJ]\>"
syn match   prmNumber   "\<\d\+[eE][+-]\=\d\+[jJ]\=\>"
syn match   prmNumber   "\<\d\+\.\%([eE][+-]\=\d\+\)\=[jJ]\=\%(\W\|$\)\@="
syn match   prmNumber   "\%(^\|\W\)\zs\d*\.\d\+\%([eE][+-]\=\d\+\)\=[jJ]\=\>"
syn region  prmString   start=+[uU]\=\z(['"]\)+ end="\z1" skip="\\\\\|\\\z1"
    \ contains=prmPythonEscape
syn match   prmSpaceError  display excludenl "\s\+$"

"---------------------------------------------------------------------------/
"  Aspect .prm builtin structures
"---------------------------------------------------------------------------/
syn keyword prmBuiltin      top bottom left right front back inner outer north west south east
syn keyword prmBuiltinBool  true false

syn match   prmSection      'subsection .*$' contains=prmSectionName
syn match   prmSection      'end\s*$'
syn match   prmSectionName  ' .*$' containedin=prmSection contained
syn match   prmVar          '^\s*set\s[^=]*' contains=prmSet
syn match   prmSet          '^\s*set' containedin=prmVar contained

"---------------------------------------------------------------------------/
"  Python/pyexpander syntax
"---------------------------------------------------------------------------/
syn keyword prmPythonStatement  containedin=prmPyexpander contained False None True as assert break continue
syn keyword prmPythonStatement  containedin=prmPyexpander contained del exec global lambda nonlocal pass print
syn keyword prmPythonStatement  containedin=prmPyexpander contained return with yield class def
syn keyword prmPythonCond       containedin=prmPyexpander contained elif else if
syn keyword prmPythonRepeat     containedin=prmPyexpander contained for while
syn keyword prmPythonOperator   containedin=prmPyexpander contained and in is not or
syn keyword prmPythonException  containedin=prmPyexpander contained except finally raise try
syn keyword prmPythonInclude    containedin=prmPyexpander contained from import
syn keyword prmPythonAsync      containedin=prmPyexpander contained async await
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained False True None
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained NotImplemented Ellipsis __debug__
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained abs all any bin bool bytearray callable chr
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained classmethod compile complex delattr dict dir
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained divmod enumerate eval filter float format
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained frozenset getattr globals hasattr hash
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained help hex id input int isinstance
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained issubclass iter len list locals map max
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained memoryview min next object oct open ord pow
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained print property range repr reversed round set
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained setattr slice sorted staticmethod str
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained sum super tuple type vars zip __import__
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained basestring cmp execfile file
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained long raw_input reduce reload unichr
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained unicode xrange
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained ascii bytes exec
syn keyword prmPythonBuiltin    containedin=prmPyexpander contained apply buffer coerce intern
syn match   prmPythonEscape     +\\[abfnrtv'"\\]+           containedin=prmString contained
syn match   prmPythonEscape     "\\\o\{1,3}"                containedin=prmString contained
syn match   prmPythonEscape     "\\x\x\{2}"                 containedin=prmString contained
syn match   prmPythonEscape     "\%(\\u\x\{4}\|\\U\x\{8}\)" containedin=prmString contained
" Python allows case-insensitive Unicode IDs: http://www.unicode.org/charts/
syn match   prmPythonEscape     "\\N{\a\+\%(\s\a\+\)*}"     containedin=prmString contained
syn match   prmPythonEscape     "\\$"                       containedin=prmString contained
syn match   prmPy               '\$py\s*'                   containedin=prmPyexpander contained
syn match   prmPyexpander       '\$..\s*(.*)'
    \ contains=ALLBUT,prmSection,prmSectionName,prmSet,prmVarName,prmBuiltinBool,prmBuiltin

"syn match   prmPyVar  '\$(\S*)'
syn match   prmPyVar  '\$([^)]*)'
syn match   prmPyCond '\$if'
syn match   prmPyCond '\$else'
syn match   prmPyCond '\$endif'

"----------------------------------------------------------------------------/
"  Setup syntax highlighting
"----------------------------------------------------------------------------/

let b:current_syntax = "prm"

" Comments
hi def link prmTodo             Todo
hi def link prmComment          Comment

" Constants
hi def link prmString           String
hi def link prmNumber           Number
hi def link prmBuiltin          Constant
hi def link prmBuiltinBool      Boolean
hi def link prmPyBuiltin        Constant
hi def link prmSpaceError       Error

" prm syntax
hi def link prmSection          Structure
hi def link prmSectionName      Function
hi def link prmSet              Operator
hi def link prmVar              Identifier

" python / pyexpander macros
hi def link prmPy               Macro
hi def link prmPyVar            Macro
hi def link prmPyCond           Macro
hi def link prmPythonStatement  Statement
hi def link prmPythonCond       Conditional
hi def link prmPythonRepeat     Repeat
hi def link prmPythonOperator   Operator
hi def link prmPythonException  Exception
hi def link prmPythonInclude    Include
hi def link prmPythonAsync      Statement
hi def link prmPythonBuiltin    Function
hi def link prmPythonEscape     Special
