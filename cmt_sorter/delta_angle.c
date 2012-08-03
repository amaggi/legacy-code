#include<math.h>

float delta_angle(float a, float b){

    return ((int)(a - b) + 360) % 360;

}

float sum_angle(float a, float b){

    return (int)(a - b) % 360;

}

