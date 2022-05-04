namespace Harmonic;
public static class Function
{
    public static uint      NumberFunc;     /// Номер задачи
    public static double    omega;          /// Значение omega

    public static void Init(uint numF) {
        NumberFunc = numF;

        switch(NumberFunc) {
            case 1: 
                omega = 1;          /// Полином первой степени
            break;
            case 2: 
                omega = 1;          /// Полином второй степени
            break;
            case 3:
                omega = 1;          /// Полином третьей степени 
            break;
            case 4:
                omega = 1;          /// Полином четвертой степени 
            break;
            case 5:
                omega = 1;          /// Неполином
            break;

            /// Исследование на большом и небольшом количестве узлов
            // case 2: 
            //     // omega = 1e-4;
            //     // omega = 1e+4;
            //     // omega = 1e+9;
            // break;
        }
    }

    //* U-функция синуса
    public static double Us(Vector vec) {
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => x + y + z,                               /// Полином первой степени
            2 => Pow(x, 2) + Pow(y, 2) + Pow(z, 2),       /// Полином второй степени
            3 => Pow(x, 3) + Pow(y, 3) + Pow(z, 3),       /// Полином третьей степени
            4 => Pow(x, 4) + Pow(y, 4) + Pow(z, 4),       /// Полином четвертой степени
            5 => Exp(-x-y-z),                             /// Неполином

            _ => 0,
        };
    }

    //* U-функция косинуса
    public static double Uc(Vector vec) {
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => x - y - z,                             /// Полином первой степени
            2 => Pow(x, 2) - Pow(y, 2) - Pow(z, 2),     /// Полином второй степени
            3 => Pow(x, 3) - Pow(y, 3) - Pow(z, 3),     /// Полином третьей степени
            4 => Pow(x, 4) - Pow(y, 4) - Pow(z, 4),     /// Полином четвертой степени
            5 => Exp(-x-y),                             /// Неполином

            _ => 0,
        };
    }

    //* F-функция от синуса
    public static double Fs(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => -omega*Sigma(vec)*Uc(vec) - omega*omega*Hi(vec)*Us(vec),                                                   /// Полином первой степени
            2 => -Lambda(vec)*6 - omega*Sigma(vec)*Uc(vec) - omega*omega*Hi(vec)*Us(vec),                                   /// Полином второй степени
            3 => -6*(x + y + z) - omega*Sigma(vec)*Uc(vec) - omega*omega*Hi(vec)*Us(vec),                                   /// Полином третьей степени
            4 => -12*(Pow(x, 2) + Pow(y, 2) + Pow(z, 2)) - omega*Sigma(vec)*Uc(vec) - omega*omega*Hi(vec)*Us(vec),          /// Полином четвертой степени
            5 => -6.01*Exp(-x-y-z) - Exp(-x-y),                                                                             /// Неполином

            _ => 0,
        };
    }  

    //* F-функция от косинуса
    public static double Fc(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => omega*Sigma(vec)*Us(vec) - omega*omega*Hi(vec)*Uc(vec),                                                     /// Полином первой степени
            2 => Lambda(vec)*2 + omega*Sigma(vec)*Us(vec) - omega*omega*Hi(vec)*Uc(vec),                                     /// Полином второй степени
            3 => -6*(x - y - z) + omega*Sigma(vec)*Us(vec) - omega*omega*Hi(vec)*Uc(vec),                                    /// Полином третьей степени
            4 => -12*(Pow(x, 2) - Pow(y, 2) - Pow(z, 2)) + omega*Sigma(vec)*Us(vec) - omega*omega*Hi(vec)*Uc(vec),           /// Полином четвертой степени
            5 => -4.01*Exp(-x-y) + Exp(-x-y-z),                                                                              /// Неполином

            _ => 0,
        };
    } 

    //* Лямбда
    public static double Lambda(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => 1,       /// Полином первой степени
            2 => 1,       /// Полином второй степени
            3 => 1,       /// Полином третьей степени
            4 => 1,       /// Полином четвертой степени
            5 => 2,       /// Неполином

            /// Исследование на большом и небольшом количестве узлов
            // 2 => 1e+2,
            // 2 => 1e+4,
            // 2 => 8e+5,

            _ => 0,
        };
    }

    //* Сигма
    public static double Sigma(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => 1,       /// Полином первой степени
            2 => 1,       /// Полином второй степени
            3 => 1,       /// Полином третьей степени
            4 => 1,       /// Полином четвертой степени
            5 => 1,       /// Неполином

            /// Исследование на большом и небольшом количестве узлов
            // 2 => 0,
            // 2 => 1e+4,
            // 2 => 1e+8,

            _ => 0,
        };
    }

    //* Кси
    public static double Hi(Vector vec) { 
        (double x, double y,  double z) = vec;
        return NumberFunc switch 
        {
            1 => 0.01,       /// Полином первой степени
            2 => 0.01,       /// Полином второй степени
            3 => 0.01,       /// Полином третьей степени
            4 => 0.01,       /// Полином четвертой степени
            5 => 0.01,       /// Неполином

            /// Исследование на большом и небольшом количестве узлов
            // 2 => 8.81e-12,
            // 2 => 1e-11,
            // 2 => 1e-10,


            _ => 0,
        };
    }



}
