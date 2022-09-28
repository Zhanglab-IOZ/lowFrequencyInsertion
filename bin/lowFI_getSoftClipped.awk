#!/usr/bin/awk -f
### 2022.09.06



BEGIN{
    marker = 1
}
{
    split($6,a,"S")
    if (NR%2 == 1){
        if (length(a) == 2){
            match($6,/([0-9]+)(S)/,a)
            if (a[1] <= up && a[1] >= low){
                printf "%s\t%s\n",$0,"XG:Z:CLIP"
                marker = 0
            }
            else {
                interline = $0
            }
        }
        else {
            interline = $0
        }
    }
    else if (marker == 0){
            if (length(a) == 2){
                match($6,/([0-9]+)(S)/,a)
                if (a[1] <= up && a[1] >= low){
                    printf "%s\t%s\n",$0,"XG:Z:CLIP"
                }
                else {
                    printf "%s\t%s\n",$0,"XG:Z:UNCLIP"
                }                
            }
            else {
                printf "%s\t%s\n",$0,"XG:Z:UNCLIP"
            }
            marker = 1
        }
    else {
        if (length(a) == 2){
            match($6,/([0-9]+)(S)/,a)
            if (a[1] <= up && a[1] >= low){
                printf "%s\t%s\n",$0,"XG:Z:CLIP"
                printf "%s\t%s\n",interline,"XG:Z:UNCLIP"
            }
        }            
    }
}