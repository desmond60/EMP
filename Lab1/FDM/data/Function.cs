namespace FDM;
public static class Function
{
    public static uint   NumberFunc;     /// Номер задачи
    public static double lambda;         /// Лямбда
    public static double gamma;          /// Гамма

    //* Инициализации лямбды и гаммы
    public static void Init() {
        switch(NumberFunc) 
        {
            case 1: /// easy
            lambda = 2; gamma = 3;
            break;

            case 2: /// sec_kraev
            lambda = 2; gamma = 3;
            break;

            case 3: /// polynom_2
            lambda = 2; gamma = 3;
            break;

            case 4: /// polynom_3
            lambda = 2; gamma = 3;
            break;

            case 5: /// polynom_4
            lambda = 2; gamma = 3;
            break;

            case 6: /// not_polynom
            lambda = 1; gamma = 1;
            break;
        }
    }

    //* Функция u(x,y)
    public static double Absolut(double x, double y, uint side = 0)
    { 
        switch(NumberFunc)
        {
            case 1: /// easy
            return 2*x + 4*y;

            case 2: /// sec_kraev
            return side switch
            {
                3 => 4,
                6 => 8,
                _ => 2*x + 4*y
            };

            case 3: /// polynom_2
            return 2*Pow(x, 2) + 4*Pow(y, 2);

            case 4: /// polynom_3
            return 2*Pow(x, 3) + 4*Pow(y, 3);

            case 5: /// polynom_4
            return 2*Pow(x, 4) + 4*Pow(y, 4);
            
            case 6: /// not_polynom
            return Sin(x + y);
        }
        return 0;
    }

    //* Функция f(x,y)
    public static double Func(double x, double y)
    { 
        switch(NumberFunc)
        {
            case 1: /// easy
            return 6*x + 12*y;

            case 2: /// sec_kraev
            return 6*x + 12*y;

            case 3: /// polynom_2
            return 6*Pow(x, 2) + 12*Pow(y, 2) - 24;

            case 4: /// polynom_3
            return -24*x - 48*y + 6*Pow(x, 3) + 12*Pow(y, 3);

            case 5: /// polynom_4
            return -48*Pow(x, 2) - 96*Pow(y, 2) + 6*Pow(x, 4) + 12*Pow(y, 4);

            case 6: /// not_polynom
            return 3*Sin(x + y);
        }
        return 0;
    }   

    //* Функция проверки имеется ли на стороне области второе краевое
    public static bool IsSecondKraev(uint side)
    {
        switch(NumberFunc)
        {
            case 1: /// easy
            return side switch
            {
                _ => false,
            };

            case 2: /// sec_kraev
            return side switch
            {
                3 => true,
                6 => true,
                _ => false
            };

            case 3: /// polynom_2
            return side switch
            {
                _ => false,
            };

            case 4: /// polynom_3
            return side switch
            {
                _ => false,
            };

            case 5: /// polynom_4
            return side switch
            {
                _ => false,
            };

            case 6: /// not_polynom
            return side switch
            {
                _ => false,
            };
        }
        return false;
    }
}
