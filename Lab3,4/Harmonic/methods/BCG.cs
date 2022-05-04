namespace Harmonic.methods;
public class BCG : LUMethods, IMethods
{
    private SLAU slau;                       /// Матрица

    private double EPS;                      /// Точность
    private int maxIter;                     /// Максимальное количество итераций 

    public string Path { get; set; }         /// Путь к задаче

    //: Табличка с итерациями и невязкой
    StringBuilder table_iter;                /// Табличка каждой итерации и невязки
    //: -----------------------------------------------------

    public BCG(SLAU slau, double EPS, int maxIter) {
        this.slau    = slau;
        this.EPS     = EPS;
        this.maxIter = maxIter;
    }

    public Vector solve(bool log = false) {
        //: Таблички
        if (log == true)
            table_iter = new StringBuilder($"Iter{" ", 5} Nev{" ", 12}\n");
        //: ---------------------------------------------------

        Vector r0   = new Vector(slau.N);
        Vector rk   = new Vector(slau.N);
        Vector z    = new Vector(slau.N);
        Vector LAUz = new Vector(slau.N);
        Vector p    = new Vector(slau.N);
        Vector LAUp = new Vector(slau.N);
        Vector y    = new Vector(slau.N);

        double f_norm = Sqrt(Scalar(slau.f, slau.f));

        // LU-разложение
        Init(slau);
        LU_decomposition();

        r0 = Solve_L(slau.f - slau.Mull(slau.q));
        Vector.Copy(r0, rk);
        Vector.Copy(slau.q, y);
        z = Solve_U(r0);

        int Iter = 0;
        double Eps = 0, alpha, gamma, betta;
        do {
            betta  = Scalar(rk, r0);
            LAUz   = Solve_L(slau.Mull(Solve_U(z)));
            alpha  = betta / Scalar(LAUz, r0);
            p      = rk - alpha*LAUz;
            LAUp   = Solve_L(slau.Mull(Solve_U(p)));
            gamma  = Scalar(p, LAUp) / Scalar(LAUp, LAUp);
            y    = y + alpha*z + gamma*p;
            rk     = p - gamma*LAUp;
            betta  = (alpha*Scalar(rk, r0)) / (gamma*betta);
            z      = rk + betta*z - betta*gamma*LAUz;
            Eps    = Sqrt(Scalar(rk, rk)) / f_norm;

            //: Табличка
            if (log == true) {
                table_iter.Append($"{String.Format("{0,4}", ++Iter)}"  + 
                                  $"{String.Format("{0,19}",  Eps.ToString("E6"))}\n");
            }
            //: -------------------------------------------------
        } while (Iter < maxIter &&
                 Eps  > EPS);
        
        //: Запись таблички в файл
        if (log == true)
            File.WriteAllText(Path + "/table_iter.txt", table_iter.ToString());
        //: -------------------------------------------------
        slau.q = Solve_U(y);
        return slau.q;
    }
    
    public string GetName() => "BCGSTAB";
}
