#!/usr/bin/awk -f
### 2022.09.06



BEGIN{
    marker = 1
}
{
    split($1,a,":")
    if (NR%2 == 1){
        if ($3 == "ALU" && a[length(a)] == "CLIP")
            interline = $0
        else {
            marker = 0
        }
    }
    else if (marker == 0){
        marker = 1
    }
    else {
        if ($3 == "ALU" && a[length(a)] == "CLIP"){
            printf "%s\n%s\n",interline,$0     
        }
    }
}
