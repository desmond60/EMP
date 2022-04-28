namespace Harmonic.methods;
public class LOS : LUMethods, IMethods
{
    private SLAU slau;                       /// Матрица
    protected double EPS;                    /// Точность решения 
    protected int maxIter;                   /// Максимальное количество итераций

    public string Path { get; set; }         /// Путь к задаче

    //: Табличка с итерациями и невязкой
    StringBuilder table_iter;                /// Табличка каждой итерации и невязки
    //: -----------------------------------------------------

    public LOS(SLAU slau, double EPS, int maxIter) {
        this.slau    = slau;
        this.EPS     = EPS;
        this.maxIter = maxIter;
    }

    public Vector solve(bool log = false) {
        //: Таблички
        if (log == true)
            table_iter = new StringBuilder($"Iter{" ", 5} Nev{" ", 12}\n");
        //: ---------------------------------------------------

        Vector r   = new Vector(slau.N);
        Vector z   = new Vector(slau.N);
        Vector Az  = new Vector(slau.N);
        Vector LAU = new Vector(slau.N);
        Vector p   = new Vector(slau.N);

        double f_norm = Sqrt(Scalar(slau.f, slau.f));

        // LU-разложение 
        Init(slau);
        LU_decomposition();

        r = Solve_L(slau.f - slau.Mull(slau.q));
        z = Solve_U(r);
        p = Solve_L(slau.Mull(z));

        double alpha, betta, Eps = 0;
        int Iter = 0;
        do {
            betta  = Scalar(p, p);
            alpha  = Scalar(p, r) / betta;
            slau.q = slau.q + alpha * z;
            r      = r - alpha * p;
            LAU    = Solve_L(slau.Mull(Solve_U(r)));
            betta  = -Scalar(p, LAU) / betta;
            z      = Solve_U(r) + betta * z;
            p      = LAU + betta * p;
            Eps    = Scalar(r, r) / f_norm;

            //: Табличка
            if (log == true) {
                table_iter.Append($"{String.Format("{0,4}", ++Iter)}"  + 
                                  $"{String.Format("{0,19}",  Eps.ToString("E6"))}\n");
            }
            //: -------------------------------------------------
        } while(Iter < maxIter &&
                Eps  > EPS);     
    

        //: Запись таблички в файл
        if (log == true)
            File.WriteAllText(Path + "/table_iter.txt", table_iter.ToString());
        //: -------------------------------------------------

        return slau.q;
    }

    public string GetName() => "LOS";
}
