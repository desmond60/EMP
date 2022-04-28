namespace Harmonic.methods;
public class LU : IMethods
{
    private SLAU slau;

    public string Path { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

    /// Матрица

    public LU(SLAU slau) {
        this.slau = slau;
    }

    public string GetName() => "LU";

    public Vector solve(bool log = false)
    {
        throw new NotImplementedException();
    }
}
