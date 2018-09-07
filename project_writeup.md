
---

**Path Planning Project**

## [Rubric](https://review.udacity.com/#!/rubrics/1020/view) Points
###Here I will consider the rubric points individually and describe how I addressed each point in my implementation.  

---
###Writeup / README

###Reflection

####1. Describe code to generate path.

This was quite an interesting start to the term and I could  talk about state machines, debugging and other interesting/challenging topics, but I'll focus on path generation.

Lines 302 to 521 in my main.cpp file includes code for determining the desired lane and speed of my car.

After this, lines 531 to 558 added initial reference/anchor points to my path to help smoothen the path to be generated.

In lines 560 to 572, I got the way points for three evenly spaced segments of 30 meters each (using my desired lane as a variable), then I added them to my path xy variables.

In lines 574 to 582, I adjusted my path points to ensure the reference/anchor angle is 0 degrees.

Finally, in lines 584 to 628, I used my desired speed and other variables to calculate the number of N points I'ld need to have my desired motion.

Using the spline library I was able to get the desired y value for each x value on my desired path.