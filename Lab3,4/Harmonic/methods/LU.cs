namespace Harmonic.methods;
public class LU : LUMethods, IMethods
{
    private SLAU slau;                       /// Матрица

    public string Path { get; set; }         /// Путь к задаче

    public LU(SLAU slau) { this.slau = slau; }

    public Vector solve(bool log = false) {
        Init_profile(slau);                                    //? СЛАУ в профильный формат из разреженного
        LU_decomposition_profile();                            //? LU-разложение
        Vector q = reverseRunning(normalRunning(slau.f));      //? Решение прямой и обратный ход
        return q;
    }

    public string GetName() => "LU";
}
