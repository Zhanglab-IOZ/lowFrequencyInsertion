#!/usr/bin/awk -f
### 2022.08.31



BEGIN{
    marker = 1
}
{
    split($6,a,"S")
    if (NR%2 == 1){
        if (length(a) == 2){
            match($6,/([0-9]+)(S)/,a)
            if (a[1] <= up && a[1] >= low){
                print $0
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
            print $0
            marker = 1
        }
    else {
        if (length(a) == 2){
            match($6,/([0-9]+)(S)/,a)
            if (a[1] <= up && a[1] >= low){
                print $0
                print interline
            }
        }            
    }
}