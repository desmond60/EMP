namespace FDM;
public class Seidel
{
    private Matrix matrix; /// Матрица 
    private double omega;  /// Параметр релаксации

    public Seidel(Matrix matrix, int iter, double eps, double omega = 1) {
        this.matrix         = matrix;
        this.omega          = omega;
        this.matrix.maxIter = iter;
        this.matrix.EPS     = eps;
    } 

    //* Решение СЛАУ
    public void solve(bool flag = false) {
        double sum, Nev = 0, norm_f;
        int Iter = 0;
        norm_f = Norm(matrix.pr);

        do {
            Nev = 0;
            for (int i = 0; i < matrix.N; i++) {
                sum = matrix.di[i] * matrix.x[i];
                
                if (i < matrix.N - 1)
                    sum += matrix.du1[i] * matrix.x[i + 1];
                if (i >= 2)
                    sum += matrix.dl1[i - 1] * matrix.x[i - 1];

                if (i < matrix.dshift)
                    sum += matrix.du2[i] * matrix.x[matrix.shift1 + i];
                else if (i < matrix.N - matrix.shift2) 
                    sum += matrix.du2[i] * matrix.x[matrix.shift2 + i];       
                
                if (i >= matrix.shift1 + matrix.dshift)
                    sum += matrix.dl2[i - matrix.shift1] * matrix.x[i - matrix.shift2];
                else if (i >= matrix.shift1)
                    sum += matrix.dl2[i - matrix.shift1] * matrix.x[i - matrix.shift1];
                
                Nev += (matrix.pr[i] - sum) * (matrix.pr[i] - sum);
                matrix.x[i] += omega / matrix.di[i] * (matrix.pr[i] - sum);
            }

            Nev = Sqrt(Nev) / norm_f; /// Относительная невязка
            Iter++;
            if (flag)
                WriteLine($"Iter: {Iter, -10} Nev: {Nev.ToString("E3")}");
        } while (Nev > matrix.EPS &&
                    Iter <= matrix.maxIter);
    }
}