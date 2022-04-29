namespace Harmonic;
public static class Function
{
    public static uint      NumberFunc;     /// Номер задачи
    public static double    omega;         

    public static void Init(uint numF) {
        NumberFunc = numF;

        switch(NumberFunc) {
            case 1: 
                omega  = 1000;
            break;
        }
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
            1 => -Hi(vec)*Us(vec)*Pow(omega, 2) - Sigma(vec)*Uc(vec)*omega - 6*Lambda(vec)*(x + y + z),
            _ => 0,
        };
    }  

    //* F-функция от косинуса
    public static double Fc(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => -Hi(vec)*Uc(vec)*Pow(omega, 2) + Sigma(vec)*Us(vec)*omega - 6*Lambda(vec)*(2*x - y + 3*z),
            _ => 0,
        };
    } 

    //* Лямбда
    public static double Lambda(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => 1000,
            _ => 0,
        };
    }

    //* Сигма
    public static double Sigma(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => 3,
            _ => 0,
        };
    }

    //* Кси
    public static double Hi(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => 1e-11,
            _ => 0,
        };
    }



}
