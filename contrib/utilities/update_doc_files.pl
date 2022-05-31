#!/bin/perl

# This script fixes some artifacts in the automatically converted documentation.
# Do not run this script directly, run it via the "update_doc_files.sh".

while (<>)
{
    # Convert embed into img
    s/\<embed/\<img/g;

    # Remove newlines in alt text
    s/(alt=.*[^>])\n/$1/gs;

    # convert img without environment into figures
    s/(?<!\<figure\>\n)\<img(.*?)\/\>/\n\<figure\>\n\<img$1\/\>\n\<\/figure\>/gs;

    # Replace .prm links
    s/\[(.*?\.prm)\]\[\]/\[$1\]\(https\:\/\/www\.github\.com\/geodynamics\/aspect\/blob\/main\/$1\)/g;

    # Replace .cc links
    s/\[(.*?\.cc)\]\[\]/\[$1\]\(https\:\/\/www\.github\.com\/geodynamics\/aspect\/blob\/main\/$1\)/g;

    # Replace directory links
    s/\[(benchmarks\/[^\s]*?)\]\[\]/\[$1\]\(https\:\/\/www\.github\.com\/geodynamics\/aspect\/blob\/main\/$1\)/g;
    s/\[(cookbooks\/[^\s]*?)\]\[\]/\[$1\]\(https\:\/\/www\.github\.com\/geodynamics\/aspect\/blob\/main\/$1\)/g;

    # Fix parameter references
    s/\[\\\[(parameters:.*)\\\]\]\[[0-9]*\]/{ref}\`$1\`/g;

    # Fix section references
    s/\[\\\[(sec:.*)\\\]\]\[[0-9]*\]/{ref}\`$1\`/g;

    # Replace sub-subsections
    s/\#\#\#\#\#\#/\#\#\#/g;

    # Replace subsections
    s/\#\#\#\#\#/\#\#/g;

    # Fix section headings
    s/\#\#\#\#/\#/g;

    # Replace an unnecessary triple backtick
    s/\&ldquo\;\`(.*?)\`\&rdquo\;/\`$1\`/g;

    # Fix ASPECT's in documentation
    s/ \&rsquo\;s/ ASPECT\'s/g;

    # Fix figure references
    s/src\=\".*\/(.*\.png)\"/src\=\"$1\"/g;
    s/src\=\".*\/(.*\.jpg)\"/src\=\"$1\"/g;
    s/src\=\".*\/(.*)\.pdf\"/src\=\"$1\.svg\"/g;
    s/src\=\".*\/(.*)\.\*\"/src\=\"$1\.\*\"/g;

    # Convert figure environments
    s/\<figure\>/\`\`\`\{figure-md\}/g;
    s/\<\/figure\>/\`\`\`/g;

    # Use alt text instead of figure caption
    s/ alt=\"(.*?)\"(.*?)\<figcaption.*?\<\/figcaption\>/$2 $1/g;
    s/ alt=\"(.*?)\"(.*?)\/\>/$2\/\> $1/g;

    # Remove figure caption tags
    s/\<figcaption aria-hidden\=\"true\"\>\<em\>//g;
    s/\<\/em\>\<\/figcaption\>//g;

    # Insert necessary newlines
    s/(\<img.*\/\>)([^\n])/$1\n\n$2/g;

    # Remove title tag in figures
    s/\stitle\=\"fig\:"//g;

    # Convert table environments
    s/\<div id\=\"(tab\:[^\s]*?)\"\>.*?(\|.*?\|)(?:\n{2})(.*?)(?:\n{2})\<\/div\>/\`\`\`\{table\} $3\n\:name\: $1\n\n$2\n\n\`\`\`/gs;
    # Remove new lines from table captions
    s/(?:\`\`\`\{table\}|(?!^)\G)(?:(?!\:name\:)[^\n]*)\K\n((?!\:name\:))/ /gs;

    # Remove empty brackets
    s/\[\]//g;

    # Replace quotes
    s/\&lsquo\;/'/g;
    s/\&rsquo\;/'/g;

    s/\&ldquo\;/"/g;
    s/\&rdquo\;/"/g;

    # Replace dashes
    s/\&ndash\;/--/g;
    s/\&mdash\;/---/g;

    # Replace invalid latex
    s/textsubscript/\_/g;

    # Extract label
    s/(\`\`\`\{figure-md\})\n(.*?)id\=\"(.*?)\" /$1 $3\n$2/gs;

    # Remove unnecessary div center around figures
    s/\<div class\=\"center\"\>\n\n(\`\`\`\{figure-md\}.*?)\<\/div\>/$1/gs;

    print;
}

