void calendar (long year, int *month, int *day, long jday)
{
int leap=0;


/* test for leap year */
if ( (year%4 == 0 && year%100 != 0) || (year%400 == 0) )
        leap=1;

if (jday < 32 ) {
        *day=jday;
        *month=1;
}
else if (jday <= 59 + leap) {
        *day = jday - 31;
        *month = 2;
}
else if (jday <= 90 + leap) {
        *day = jday - (59 + leap);
        *month = 3;
}
else if (jday <= 120 + leap) {
        *day = jday - (90+leap);
        *month=4;
}
else if (jday <= 151 + leap) {
        *day = jday - (120 + leap);
        *month = 5;
}
else if (jday <= 181 + leap) {
        *day = jday - (151 + leap);
        *month = 6;
}else if (jday <= 212 + leap) {
        *day = jday - (181 + leap);
        *month = 7;
}
else if (jday <= 243 + leap) {
        *day = jday - (212 + leap);
        *month = 8;
}
else if (jday <= 273 + leap) {
        *day = jday - (243 + leap);
        *month = 9;
}
else if (jday <= 304 + leap) {
        *day = jday - (273 + leap);
        *month = 10;
}
else if (jday <= 334 + leap) {
        *day = jday - (304 + leap);
        *month = 11;
}
else {
        *day = jday - (334 + leap);
        *month = 6;
}

}
