namespace NonLinear;
public static class Function
{
    public static uint   NumberFunc;     /// Номер задачи

    //* Функция u(x,y)
    public static double Absolut(double x)
    { 
        switch(NumberFunc)
        {
            case 1: /// test1-evenly (const)
            return 2.5;

            case 2: /// test2-polynom1
            return x;

            case 3: /// test3-polynom2
            return x*x;

            case 4: /// test4-polynom3
            return Pow(x, 3);

            case 5: /// test5_nopolynom 
            return Sin(2*x);
        }
        return 0;
    }

    //* Функция f(x,y)
    public static double Func(double x)
    { 
        switch(NumberFunc)
        {
            case 1: /// test1-evenly (const)
            return 5;

            case 2: /// test2-polynom1
            return 2*x;

            case 3: /// test3-polynom2
            return 2*x*x -8*x -8;

            case 4: /// test4-polynom3
            return -34*Pow(x, 3) - 24*x;

            case 5: /// test5_nopolynom 
            return 16*Sin(2*x)*Cos(2*x);
        }
        return 0;
    }  

    //* Лямбда от производной l(u`)
    public static double Lambda(double diff, double x)
    { 
        switch(NumberFunc)
        {
            case 1: /// test1-evenly (const)
            return diff + 4;

            case 2: /// test2-polynom1
            return diff + 4;

            case 3: /// test3-polynom2
            return diff + 4;

            case 4: /// test4-polynom3
            return diff + 4;

            case 5: /// test5_nopolynom
            return diff + 4;
        }
        return 0;
    } 

    //* Производная от Лямбды (в зависимости почему берем)
    public static double DiffLambda(double q, double q_next, double diff, double x, double h,  uint num_q)
    { 
        switch(NumberFunc)
        {
            case 1: /// test1-evenly (const)
            return num_q switch 
            {
                1 => -1/h,
                2 =>  1/h,
                _ => 0
            };

            case 2: /// test2-polynom1
            return num_q switch 
            {
                1 => -1/h,
                2 =>  1/h,
                _ => 0
            };

            case 3: /// test3-polynom2
            return num_q switch 
            {
                1 => -1/h,
                2 =>  1/h,
                _ => 0
            };

            case 4: /// test4-polynom3
            return num_q switch 
            {
                1 => -1/h,
                2 =>  1/h,
                _ => 0
            };

            case 5: /// test5_nopolynom
            return num_q switch 
            {
                1 => -1/h,
                2 =>  1/h,
                _ => 0
            };
        }
        return 0;
    }

    //* Производная абсолютной функции
    public static double Diff(double x)
    { 
        switch(NumberFunc)
        {
            case 1: /// test1-evenly (const)
            return 0;

            case 2: /// test2-polynom1
            return 1;

            case 3: /// test3-polynom2
            return 2*x;

            case 4: /// test4-polynom3
            return 3*Pow(x, 2);

            case 5: /// test5_nopolynom
            return Cos(2*x);
        }
        return 0;
    } 

}
