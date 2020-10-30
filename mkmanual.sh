#!/bin/sh

cat |
    sed -n '{/NAME/,/AUTHORS/p;/AUTHORS/{n;p}}' |
    sed -E 's/Ss\{([^}]*)\}/## \1\n/g' |
    sed -E 's/Nm\{([^}]*)\}/`\1`/g' |
    sed -E 's/Ql\{([^}]*)\}/`\1`/g' |
    sed -E 's/Ql\[([^]]*)\]/`\1`/g' |
    sed -E 's/ *Dl\{([^}]*)\}/>***\1***/g' |
    sed -E 's/Cm\{([^}]*)\}/**\1**/g' |
    sed -E 's/Fl\{([^}]*)\}/\1/g' |
    sed -E 's/Ar\{([^}]*)\}/*\1*/g' |
    sed -E 's/Ev\{([^}]*)\}/`\1`/g' |
    sed -E 's/^    //' |
    sed -E 's/\\$/\\\\\\/' |
    sed -E '/COMMANDS/,/##/{s/^([^ #].*)/* \1\n/;s/^    /  /}' |
    sed -E '/OPTIONS/,/##/{s/^([^ #].*)  (.*)/* \1\n\n  \2\n/;s/^    /  /}' |
    sed -E '/ARGUMENTS/,/##/{s/^([^ #].*)  (.*)/* \1\n\n  \2\n/;s/^    /  /}' |
    sed -E '/AUTHORS/,$s/<[^>]+>//g' |
    sed -E 's/ +$//g' |
    cat -s
