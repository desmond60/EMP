namespace Harmonic;
public static class Function
{
    public static uint      NumberFunc;     /// Номер задачи
    public static double    omega;
    public static double    hi;
    public static double    sigma; 
    public static double    lambda;          

    public static void Init(uint numF) {
        NumberFunc = numF;

        switch(NumberFunc) {
            case 1: 
                omega  = 1000;
                lambda = 1000;
                sigma  = 3;
                hi     = 1e-11;
            break;
        }
    }

    //* Абсолютное значение функции U
    public static double Absolut(double x) { 
        switch(NumberFunc)
        {
            case 1: /// test1-evenly (const)
            return 2.5;

        }
        return 0;
    }

    //* U-функция синуса
    public static double Us(Vector vec) {
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => Pow(x, 3) + Pow(y, 3) + Pow(z, 3),
            _ => 0,
        };
    }

    //* U-функция косинуса
    public static double Uc(Vector vec) {
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => 2*Pow(x, 3) - Pow(y, 3) + 3*Pow(z, 3),
            _ => 0,
        };
    }

    //* F-функция от синуса
    public static double Fs(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => -hi*Us(vec)*Pow(omega, 2) - sigma*Uc(vec)*omega - 6*lambda*(x + y + z),
            _ => 0,
        };
    }  

    //* F-функция от косинуса
    public static double Fc(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => -hi*Uc(vec)*Pow(omega, 2) - sigma*Uc(vec)*omega - 6*lambda*(2*x - y + 3*z),
            _ => 0,
        };
    }  


}
