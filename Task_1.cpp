#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

int algoritm = 1;

class ErrorNoFile
{
    string str = "ErrorNoFile";

public:
    ErrorNoFile(string s) : str(s) {}
    void Message()
    {
        cout << "ErrorNoFile " << str << endl;
    }
};

class ErrorInterpolation
{
    string msg;

public:
    ErrorInterpolation(string s) : msg(s) {}
    void Message()
    {
        cout << "ErrorInterpolation: " << msg << endl;
    }
};

class ErrorDivideByZero
{
    string msg;

public:
    ErrorDivideByZero(string s) : msg(s) {}
    void Message()
    {
        cout << "ErrorDivideByZero: " << msg << endl;
    }
};

double T(double x)
{
    ifstream is;
    if (abs(x) <= 1)
        is.open("dat_X_1_1.dat");
    else if (x < -1)
    {
        x = 1 / x;
        is.open("dat_X00_1.dat");
    }
    else if (x > 1)
    {
        x = 1 / x;
        is.open("dat_X1_00.dat");
    }

    if (!is)
        throw ErrorNoFile("File not open");

    double xi, ti, xi1, ti1, dummy, t = 0;
    bool found = false;

    if (!(is >> xi1 >> ti1 >> dummy))
    {
        throw ErrorNoFile("Invalid file format or empty file");
    }

    if (xi1 == x)
        t = ti1;
    else
    {
        while (!is.eof())
        {
            xi = xi1;
            ti = ti1;
            is >> xi1 >> ti1 >> dummy;

            if ((xi < x && x < xi1) || (xi1 < x && x < xi))
            {
                if (xi1 == xi)
                    throw ErrorDivideByZero("Duplicate x-values in interpolation at x = " + to_string(xi));
                t = ti + (ti1 - ti) * (x - xi) / (xi1 - xi);
                found = true;
                break;
            }
            if (xi1 == x)
            {
                t = ti1;
                break;
            }
        }

        if (!found && xi1 != x)
            throw ErrorInterpolation("No matching interval for x = " + to_string(x));
    }
    is.close();
    return t;
}

double U(double x)
{
    ifstream is;
    if (abs(x) <= 1)
        is.open("dat_X_1_1.dat");
    else if (x < -1)
    {
        x = 1 / x;
        is.open("dat_X00_1.dat");
    }
    else if (x > 1)
    {
        x = 1 / x;
        is.open("dat_X1_00.dat");
    }

    if (!is)
        throw ErrorNoFile("File not open");
    double xi, ui, xi1, ui1, dummy, u = 0;
    bool found = false;

    if (!(is >> xi1 >> dummy >> ui1))
    {
        throw ErrorNoFile("Invalid file format or empty file");
    }

    if (xi1 == x)
        u = ui1;
    else
    {
        while (!is.eof())
        {
            xi = xi1;
            ui = ui1;
            is >> xi1 >> dummy >> ui1;

            if ((xi < x && x < xi1) || (xi1 < x && x < xi))
            {
                if (xi1 == xi)
                    throw ErrorDivideByZero("Duplicate x-values in interpolation at x = " + to_string(xi));
                u = ui + (ui1 - ui) * (x - xi) / (xi1 - xi);
                found = true;
                break;
            }
            if (xi1 == x)
            {
                u = ui1;
                break;
            }
        }

        if (!found && xi1 != x)
            throw ErrorInterpolation("No matching interval for x = " + to_string(x));
    }
    is.close();
    return u;
}

double Srz(double x, double y, double z)
{
    if (x > y)
        return T(x) + U(z) - T(y);
    else
        return T(y) + U(y) - U(z);
}

// Алгоритм 2
double Srs1(double x, double y, double z)
{
    if (z > y)
        return Srz(x, y, z) + 1.44 * y * z;
    else
        return y + 1.44 * Srz(z, x, y);
}

double Qrz1(double x, double y)
{
    if (abs(y) < 1)
        return x * Srs1(x, y, x);
    else
        return y * Srs1(y, x, y);
}

// Алгоритм 3
double Srs2(double x, double y, double z)
{
    if (z > y)
        return Srz(x, y, z) + y * x;
    else
        return y * z + Srz(z, x, y);
}

double Qrz2(double x, double y)
{
    if (abs(x) < 1)
        return x * Srs2(x, y, x);
    else
        return y * Srs2(y, x, y);
}

double Srs(double x, double y, double z)
{
    if (z > y)
    {
        if (z * z + x * y > 0)
            return Srz(x, y, z) + y * sqrt(z * z + x * y);
        else
        {
            algoritm = 2;
            return Srs1(x, y, z);
        }
    }
    else
    {
        if (x * x + z * y > 0)
            return y + Srz(z, x, y) * sqrt(x * x + z * y);
        else
        {
            algoritm = 3;
            return Srs2(x, y, z);
        }
    }
}

double Qrz(double x, double y)
{
    if (abs(x) < 1)
        return x * Srs(x, y, x);
    else
        return y * Srs1(y, x, y);
}

double Rrz(double x, double y, double z)
{
    if (algoritm == 1)
    {
        if (x > y)
            return x * z * Qrz(y, z);
        else
            return y * x * Qrz(x, y);
    }
    else if (algoritm == 2)
    {
        if (x > y)
            return x * y * Qrz1(y, z);
        else
            return x * z * Qrz1(x, y);
    }
    else
    {
        if (x > y)
            return x * y * Qrz2(y, z);
        else
            return y * z * Qrz2(x, y);
    }
}

double Grs(double x, double y, double z)
{
    return 0.1389 * Rrz(x, y, y) + 1.8389 * Rrz(x - y, z, y);
}

double fun(double x, double y, double z)
{
    return x * Grs(x, y, z) + y * Grs(x, z, y);
}

int main()
{
    double x, y, z, f;
    /*
        ofstream  outdat("dat_X_1_1.dat");
        outdat << -1.000 << "\t" << -4.935 << "\t" << 1.935 << endl;
        outdat << -0.900 << "\t" << -3.013 << "\t" << 0.464 << endl;
        outdat << -0.800 << "\t" << -2.316 << "\t" << 1.327 << endl;
        outdat << -0.700 << "\t" << -1.819 << "\t" << 1.976 << endl;
        outdat << -0.600 << "\t" << -1.425 << "\t" << 2.502 << endl;
        outdat << -0.500 << "\t" << -1.097 << "\t" << 2.951 << endl;
        outdat << -0.400 << "\t" << -0.816 << "\t" << 3.344 << endl;
        outdat << -0.300 << "\t" << -0.571 << "\t" << 3.695 << endl;
        outdat << -0.200 << "\t" << -0.357 << "\t" << 4.013 << endl;
        outdat << -0.100 << "\t" << -0.167 << "\t" << 4.303 << endl;
        outdat <<  0.000 << "\t" <<  0.000 << "\t" << 4.571 << endl;
        outdat <<  0.100 << "\t" <<  0.147 << "\t" << 4.618 << endl;
        outdat <<  0.200 << "\t" <<  0.276 << "\t" << 4.645 << endl;
        outdat <<  0.300 << "\t" <<  0.386 << "\t" << 4.652 << endl;
        outdat <<  0.400 << "\t" <<  0.477 << "\t" << 4.636 << endl;
        outdat <<  0.500 << "\t" <<  0.548 << "\t" << 4.596 << endl;
        outdat <<  0.600 << "\t" <<  0.597 << "\t" << 4.524 << endl;
        outdat <<  0.700 << "\t" <<  0.617 << "\t" << 4.412 << endl;
        outdat <<  0.800 << "\t" <<  0.597 << "\t" << 4.240 << endl;
        outdat <<  0.900 << "\t" <<  0.505 << "\t" << 3.956 << endl;
        outdat <<  1.000 << "\t" <<  0.000 << "\t" << 3.000 << endl;
        outdat.close();
    */

    /*
        ofstream outdat("dat_X1_00.dat");
        outdat << 0.000 << "\t" << -4.935 << "\t" << 1.935 << endl;
        outdat << 0.050 << "\t" << -2.663 << "\t" << 1.885 << endl;
        outdat << 0.100 << "\t" << -1.618 << "\t" << 1.834 << endl;
        outdat << 0.150 << "\t" << -0.773 << "\t" << 1.784 << endl;
        outdat << 0.200 << "\t" << -0.034 << "\t" << 1.732 << endl;
        outdat << 0.250 << "\t" <<  0.635 << "\t" << 1.679 << endl;
        outdat << 0.300 << "\t" <<  1.253 << "\t" << 1.625 << endl;
        outdat << 0.350 << "\t" <<  1.829 << "\t" << 1.570 << endl;
        outdat << 0.400 << "\t" <<  2.369 << "\t" << 1.512 << endl;
        outdat << 0.450 << "\t" <<  2.877 << "\t" << 1.452 << endl;
        outdat << 0.500 << "\t" <<  3.356 << "\t" << 1.388 << endl;
        outdat << 0.550 << "\t" <<  3.806 << "\t" << 1.322 << endl;
        outdat << 0.600 << "\t" <<  4.228 << "\t" << 1.251 << endl;
        outdat << 0.650 << "\t" <<  4.622 << "\t" << 1.175 << endl;
        outdat << 0.700 << "\t" <<  4.987 << "\t" << 1.093 << endl;
        outdat << 0.750 << "\t" <<  5.320 << "\t" << 1.003 << endl;
        outdat << 0.800 << "\t" <<  5.618 << "\t" << 0.905 << endl;
        outdat << 0.850 << "\t" <<  5.876 << "\t" << 0.796 << endl;
        outdat << 0.900 << "\t" <<  6.080 << "\t" << 0.675 << endl;
        outdat << 0.950 << "\t" <<  6.199 << "\t" << 0.536 << endl;
        outdat << 1.000 << "\t" <<  5.890 << "\t" << 0.377 << endl;
        outdat.close();
    */

    /*
        ofstream outdat("dat_X00_1.dat");
        outdat <<  0.000 << "\t" << -4.935 << "\t" << 1.935 << endl;
        outdat << -0.050 << "\t" << -4.435 << "\t" << 1.835 << endl;
        outdat << -0.100 << "\t" << -3.936 << "\t" << 1.735 << endl;
        outdat << -0.150 << "\t" << -3.440 << "\t" << 1.636 << endl;
        outdat << -0.200 << "\t" << -2.948 << "\t" << 1.537 << endl;
        outdat << -0.250 << "\t" << -2.461 << "\t" << 1.440 << endl;
        outdat << -0.300 << "\t" << -1.980 << "\t" << 1.344 << endl;
        outdat << -0.350 << "\t" << -1.506 << "\t" << 1.249 << endl;
        outdat << -0.400 << "\t" << -1.041 << "\t" << 1.156 << endl;
        outdat << -0.450 << "\t" << -0.585 << "\t" << 1.065 << endl;
        outdat << -0.500 << "\t" << -0.141 << "\t" << 0.976 << endl;
        outdat << -0.550 << "\t" <<  0.292 << "\t" << 0.889 << endl;
        outdat << -0.600 << "\t" <<  0.712 << "\t" << 0.806 << endl;
        outdat << -0.650 << "\t" <<  1.117 << "\t" << 0.724 << endl;
        outdat << -0.700 << "\t" <<  1.507 << "\t" << 0.646 << endl;
        outdat << -0.750 << "\t" <<  1.882 << "\t" << 0.572 << endl;
        outdat << -0.800 << "\t" <<  2.239 << "\t" << 0.500 << endl;
        outdat << -0.850 << "\t" <<  2.578 << "\t" << 0.432 << endl;
        outdat << -0.900 << "\t" <<  2.898 << "\t" << 0.368 << endl;
        outdat << -0.950 << "\t" <<  3.199 << "\t" << 0.308 << endl;
        outdat << -1.000 << "\t" <<  3.480 << "\t" << 0.252 << endl;
        outdat.close();
    */

    cout << " Input x y z: ";
    cin >> x >> y >> z;
    try
    {
        f = fun(x, y, z);
    }
    catch (ErrorNoFile &e)
    {
        e.Message();
        f = 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;
    }
    catch (ErrorInterpolation &e)
    {
        e.Message();
        f = 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;
    }
    catch (ErrorDivideByZero &e)
    {
        e.Message();
        f = 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;
    }
    cout << "\n fun = " << f << endl;
    return 0;
}