#!/usr/bin/awk -f

BEGIN{
    subsum = 0;
    subsum_line = 0;
    n_subsum = 0;
    grandsum = 0;
    OFS = "\t";
}
{
    if (/^0$/) {
        if (subsum > 0) {
            n_subsum++;
            print "subsum", n_subsum, subsum_line, NR - 1, subsum;
            subsum = 0;
            subsum_line = 0;
        }
    } else {
        subsum += $1;
        if (subsum_line == 0)
            subsum_line = NR;
        grandsum += $1;
    }
}
END{
    if (subsum > 0) {
        print "subsum", n_subsum, subsum_line, NR - 1, subsum;
        n_subsum++;
    }
    print "grandsum", n_subsum, 0, 0, grandsum;
}
