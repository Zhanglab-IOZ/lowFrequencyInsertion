#!/usr/bin/awk -f

function abs(num) {
    return sqrt(num*num)
}
BEGIN{
    sum_locPass = 0
    sum_locChange = 0
}
{
    if ($NF == "RA:Z:Y"){
        sum_locChange = sum_locChange + 1
    }
    if ($3!="ALU"){
        if (NR%4==2){
            inter_chr = $3
            inter_loc = $4
        }
        else if(NR%4 == 0 && inter_chr == $3 && abs(inter_loc - $4)<=600){
            sum_locPass = sum_locPass + 1
        }
    }
}
END{
    # printf "%s%d%s%s\n","Number of read pairs pass the cutoff after realignment: ",sum_locPass,", proportion: ",sum_locPass/(NR/2)
    printf "%s%d%s%s\n","Number of read realigned to a different location: ",sum_locChange,", proportion: ",sum_locChange/NR
}
