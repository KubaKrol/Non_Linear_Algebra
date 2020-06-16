#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cmath>

#define e 1.e-8

//funkcja sinus
double sinf(double x) {
    return sin(x / 4.) * sin(x / 4.) - x;
}

//pierwsza pochodna funkcji sinus
double pochodnasin(double x) {
    return 0.5 * sin(x / 4.) * cos(x / 4.) - 1.;
}

//funkcja tangens
double tanf(double x) {
    return tan(2.0 * x) - x - 1.0;
}

//pierwsza pochodna funkcji tangens
double pochodnatan(double x) {
    return 2. / (cos(2. * x) * cos(2. * x)) - 1.;
}

//Rownanie sinusa przeksztalcone do postaci x = fi(x)
double fisin(double x) {
    return sin(x / 4.) * sin(x / 4.);
}

//Rownanie tangensa przeksztalcone do postaci x = fi(x)
double fitan(double x) {
    return tan(2. * x) - 1.;
}

double bisekcja(double (*f)(double), double a, double b, int n)
{
    double srodek = 0;
    double residuum = 0;
    double estymator = 0;

    for (int i = 0; i < n; i++) 
    {
        if (f(a) * f(b) > 0) //Kontrola roznicy znakow
        {
            std::cout << "Brak miejsca zerowego w podanym przedziale" << std::endl;
            return -1;
        }
        else if (f(a) * f(b) < 0) 
        {
            srodek = (a + b) / 2;

            estymator = fabs((b - a) / 2);
            residuum = fabs(f(srodek));

            //Wybor przedzialu o przeciwnym znaku
            if ((f(a) * f(srodek)) < 0) 
            {
                b = srodek;
            }
            else 
            {
                a = srodek;
            }
        }

        std::cout << "iteracja " << i + 1 << " x: " << srodek;
        std::cout << " residuum: " << residuum << " estymator: " << estymator << std::endl;

        if (estymator <= e || residuum <= e)
        {
            return srodek;
        }
    }

    return srodek;
}

double newton(double (*f)(double), double (*fp)(double), double x, int n)
{
    double estymator = 0;
    double residuum = 0;
    double tempX;

    if (fp(x) == 0)
    {
        std::cout << "dzielenie przez 0" << std::endl;
        return -1;
    }

    for (int i = 0; i < n; i++)
    {
        double xn = x;
        double func = f(xn);
        x = xn - func / fp(xn); //Miejsce zerowe stycznej

        estymator = fabs(xn - x);
        residuum = fabs(func);

        std::cout << "iteracja " << i + 1 << " x: " << xn;
        std::cout << " residuum: " << residuum << " estymator: " << estymator << std::endl;

        if (residuum < e || estymator < e)
        {
            return xn;
        }
    }

    return 0;
}

//Metoda iteracyjna dwuargumentowa
//Zaleta - f(x) nie musi by? ró?niczkowalne
double sieczne(double (*f)(double), double xn, double xn1, int n) 
{
    double estymator;
    double residuum;
    double xn2 = 0;

    for (int i = 0; i < n; i++) 
    {
        xn2 = xn1 - f(xn1) / ((f(xn1) - f(xn)) / (xn1 - xn));
        xn1 = xn;
        xn = xn2;
        estymator = fabs(xn - xn1);
        residuum = fabs(f(xn2));

        std::cout << "iteracja " << i + 1 << " x: " << xn2;
        std::cout << " residuum: " << residuum << " estymator: " << estymator << std::endl;

        if (residuum <= e || estymator <= e) 
        {
            return xn2;
        }
    }

    return xn2;
}

double picard(double (*fi)(double), double (*f)(double), double x0, int n)
{
    double estymator = 0;
    double residuum = 0;
    double x = x0;

    for (int i = 0; i < n; i++) 
    {
        x0 = x;
        x = fi(x);

        residuum = fabs(f(x));
        estymator = fabs(x - x0);

        std::cout << "iteracja " << i + 1 << " x: " << x;
        std::cout << " residuum: " << residuum << " estymator: " << estymator << std::endl;

        if (estymator <= e || residuum < e)
        {
            return x;
        }
    }

    return x;
}

int main() 
{
    double a, b, x0, f0;

    std::cout << "bisekcja sin: " << std::endl;
    std::cout << bisekcja(sinf, -2, 1., 30) << std::endl;

    std::cout << "bisekcja tan: " << std::endl;
    std::cout << bisekcja(tanf, -0.5, 0.5, 30) << std::endl;

    std::cout << "sieczne sin: " << std::endl;
    std::cout << sieczne(sinf, 2, 3, 10) << std::endl;

    std::cout << "sieczne tan: " << std::endl;
    std::cout << sieczne(tanf, 0.5, 0.6, 10) << std::endl;

    std::cout << "picard sin: " << std::endl;
    std::cout << picard(fisin, sinf, 2., 10) << std::endl;

    std::cout << "newton sin: " << std::endl;
    std::cout << newton(sinf, pochodnasin, -1, 10) << std::endl;

    std::cout << "newton tan: " << std::endl;
    std::cout << newton(tanf, pochodnatan, 0.5, 10) << std::endl;
    return 0;
}
