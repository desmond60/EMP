namespace NonLinear.methods;
public class Newton : IMethods
{
    double w;                 /// Параметр релаксации
    private SLAU    slau;     /// структура СЛАУ
    private Data    data;     /// Данные задачи
    private Elem[] elems;     /// Конечные элементы

    public string GetName() => "Newton";

    public Newton(SLAU _slau, Data _data, Elem[] _elems, double _w = 1) {
        slau         = _slau;
        data         = _data;
        elems        = _elems;
        w            = _w;
        slau.MaxIter = 1000;
        slau.EPS     = 1e-10;
    }

    private bool Check(int iter, bool log = true) {        
        global();
        slau.Mult(slau.q, out double[] vec_nev);
        vec_nev = slau.f.Zip(vec_nev, (f, s) => f - s).ToArray();
        double value = Norm(vec_nev) / Norm(slau.f);

        if (log) 
            WriteLine($"{iter,5} : {value,10}");

        if (iter  > slau.MaxIter) return false;
        if (value < slau.EPS)     return false;

        return true;
    }

    void AddNewton(ref double[] mat_newton, int i) {
        double lam_arg   = (slau.q[i + 1] - slau.q[i]) / elems[i].h;
        double diff_l_q1_x1 = DiffLambda(slau.q[i], slau.q[i + 1], lam_arg, elems[i].x1, elems[i].h, 1);
        double diff_l_q1_x2 = DiffLambda(slau.q[i], slau.q[i + 1], lam_arg, elems[i].x2, elems[i].h, 1);
        double diff_l_q2_x1 = DiffLambda(slau.q[i], slau.q[i + 1], lam_arg, elems[i].x1, elems[i].h, 2);
        double diff_l_q2_x2 = DiffLambda(slau.q[i], slau.q[i + 1], lam_arg, elems[i].x2, elems[i].h, 2);

        // Матрица 1x4
        // A11_q1=A22_q1
        // A11_q2=A22_q2
        // A12_q1=A21_q1
        // A12_q2=A21_q2
        mat_newton[0] =  (diff_l_q1_x1 + diff_l_q1_x2) / (2 * elems[i].h); 
        mat_newton[1] =  (diff_l_q2_x1 + diff_l_q2_x2) / (2 * elems[i].h);
        mat_newton[2] = -(diff_l_q1_x1 + diff_l_q1_x2) / (2 * elems[i].h);
        mat_newton[3] = -(diff_l_q2_x1 + diff_l_q2_x2) / (2 * elems[i].h);
    }

    //: ***** Реализация функций интерфейса ***** :\\
    public string solve() {
        int Iter = 0; //? Количество итераций
        
        do {
            Array.Copy(slau.q, slau.q1, slau.q.Length);     //? q1 = q
            global();                                       //? Генерация глобальной матрицы    
            slau.q = slau.Progonka();                       //? Решаем СЛАУ методом Прогонки
            slau.q = newApprox();                           //? Получаем новое приближение (релаксация)
        } while (Check(++Iter));

        return output().ToString();
    }

    public void global() {
        slau.Clear();                                            // Очищаем СЛАУ
        for (int index = 0; index < data.countEl; index++) {     // Проход по всем КЭ
            (double[][] mat, double[] vec) = local(index);       // Генерируем локальную матрицу и вектор
            add_to_global(mat, vec, index);                      // Заносим локальное в глобальное 
        }
        kraev();                                                 // Учет краевых
    }

    public (double[][], double[]) local(int i) {
        var local_vec  = new double[2];
        var local_mat  = new double[2][];
        var mat_newton = new double[4]; 
        for (int j = 0; j < 2; j++) local_mat[j] = new double[2];

        double lam_arg = (slau.q[i + 1] - slau.q[i]) / elems[i].h;
        AddNewton(ref mat_newton, i);

        // + Матрица жесткости
        double value = (Lambda(lam_arg, elems[i].x1) + Lambda(lam_arg, elems[i].x2)) / (2*elems[i].h);
        local_mat[0][0] =  value + (mat_newton[0]*slau.q[i] + mat_newton[2]*slau.q[i + 1]);
        local_mat[0][1] = -value + (mat_newton[1]*slau.q[i] + mat_newton[3]*slau.q[i + 1]);
        local_mat[1][0] = -value + (mat_newton[2]*slau.q[i] + mat_newton[0]*slau.q[i + 1]);
        local_mat[1][1] =  value + (mat_newton[3]*slau.q[i] + mat_newton[1]*slau.q[i + 1]);

        // + Матрица масс
        value = (data.gamma * elems[i].h) / 6.0;
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                local_mat[j][k] += j == k 
                                ?  2*value
                                :    value;

        // Правая часть
        var f = new double[2] { Func(elems[i].x1), Func(elems[i].x2) };
        double add_newton_b0 = slau.q[i]     * (mat_newton[0]*slau.q[i] + mat_newton[1]*slau.q[i + 1]) + 
                               slau.q[i + 1] * (mat_newton[2]*slau.q[i] + mat_newton[3]*slau.q[i + 1]);
        double add_newton_b1 = slau.q[i]     * (mat_newton[2]*slau.q[i] + mat_newton[3]*slau.q[i + 1]) + 
                               slau.q[i + 1] * (mat_newton[0]*slau.q[i] + mat_newton[1]*slau.q[i + 1]);
        
        local_vec[0] = (elems[i].h / 6.0) * (2*f[0] + f[1]) + add_newton_b0;
        local_vec[1] = (elems[i].h / 6.0) * (2*f[1] + f[0]) + add_newton_b1;



        return (local_mat, local_vec);
    }

    public void add_to_global(double[][] mat, double[] vec, int i) {
        slau.di[i]     += mat[0][0];
        slau.di[i + 1] += mat[1][1];
        slau.du[i]     += mat[0][1];
        slau.dl[i + 1] += mat[1][0];
        slau.f[i]      += vec[0]; 
        slau.f[i + 1]  += vec[1];
    }

    public void kraev() {
        for (int i = 0; i < data.kraev.Length; i++) {
            switch (data.kraev[i])
            {
                case 1: firstKraev(i);  break;
                case 2: secondKraev(i); break;
                case 3: thirdKraev(i);  break;
            }
        }
    }

    public void firstKraev(int i) {
        switch(i) {
            // Левая граница
            case 0: 
                slau.di[0] = 1;
                slau.du[0] = 0;
                slau.f[0]  = Absolut(elems[0].x1);
            break;

            // Правая граница
            case 1:
                slau.di[data.countEl] = 1;
                slau.dl[data.countEl] = 0;
                slau.f [data.countEl] = Absolut(elems[data.countEl - 1].x2);
            break;
        }
    }

    public void secondKraev(int i) {
        switch(i) {
            // Левая граница
            case 0: 
                double lam_arg = (slau.q[1] - slau.q[0]) / elems[0].h;
                slau.f[0] -= Lambda(lam_arg, elems[0].x1) * Diff(elems[0].x1);
            break;

            // Правая граница
            case 1:
                lam_arg = (slau.q[data.countEl] - slau.q[data.countEl - 1]) / elems[data.countEl - 1].h;
                slau.f[data.countEl] += Lambda(lam_arg, elems[data.countEl - 1].x2) * Diff(elems[data.countEl - 1].x2);
            break;
        }
    }

    public void thirdKraev(int i)
    {
        switch(i) {
            // Левая граница
            case 0:
                slau.di[0] += data.betta;
                double lam_arg = (slau.q[1] - slau.q[0]) / elems[0].h;
                double res = Lambda(lam_arg, elems[0].x1) * (-1) * Diff(elems[0].x1);
                slau.f[0] += data.betta * (res + data.betta * Absolut(elems[0].x1));
            break;

            // Правая граница
            case 1:
                slau.di[data.countEl] += data.betta;
                lam_arg = (slau.q[data.countEl] - slau.q[data.countEl - 1]) / elems[data.countEl - 1].h;
                res = Lambda(lam_arg, elems[data.countEl - 1].x2) * Diff(elems[data.countEl - 1].x2);
                slau.f[data.countEl] += data.betta * (res + data.betta * Absolut(elems[data.countEl - 1].x2));
            break;
        }      
    }

    public double[] newApprox() {
        return slau.q.Zip(slau.q1, (f, s) => w*f + (1 - w)*s).ToArray();
    }

    public StringBuilder output() {
        var table = new StringBuilder();
        
        string margin = String.Join("", Enumerable.Repeat("-", 16));

        table.Append(String.Join("", Enumerable.Repeat("-", 52)) + "\n");
        table.Append($"|U{" ",-14} | U`{" ",-12} | |U`- U| {" ",-7}|\n");
        table.Append($"|" + margin + "|" + margin + "|" + margin + "|\n");

        for (int i = 0; i < data.countEl; i++)
            table.Append($"|{String.Format("{0,16}", slau.q[i].ToString("E6"))}" + 
                        $"|{String.Format("{0,16}", Absolut(elems[i].x1).ToString("E6"))}"   + 
                        $"|{String.Format("{0,16}", Abs(Absolut(elems[i].x1) - slau.q[i]).ToString("E6"))}|\n");

        table.Append($"|{String.Format("{0,16}", slau.q[data.countEl].ToString("E6"))}" + 
                     $"|{String.Format("{0,16}", Absolut(elems[data.countEl - 1].x2).ToString("E6"))}"   + 
                     $"|{String.Format("{0,16}", Abs(Absolut(elems[data.countEl - 1].x2) - slau.q[data.countEl]).ToString("E6"))}|\n");

        table.Append(String.Join("", Enumerable.Repeat("-", 52)) + "\n");
        
        return table;
    }
}
