#!/usr/bin/awk -f
### 2022.09.13



{
    printf "%s%s\n","@",$1
    split($6,a,"S")
    match($6,/([0-9]+)(S)/,b)
    if(a[2]) {
        len = b[1]
        printf "%s\n%s\n%s\n",substr($10,1,len),"+",substr($11,1,len)
    }
    else {
        len = b[1]
        printf "%s\n%s\n%s\n",substr($10,length($10)-len+1,len),"+",substr($11,length($11)-len+1,len)
    }

}
