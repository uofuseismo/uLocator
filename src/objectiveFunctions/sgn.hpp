#ifndef SGN_HPP
#define SGN_HPP
namespace
{
double sgn(const double x)
{
    if (x == 0){return 0;} 
    if (x > 0)
    {   
        return 1;
    }   
    else //if (x < 0)
    {   
        return -1; 
    }   
}
}
#endif
