namespace Harmonic.methods;
public interface IMethods
{
    public string Path { get; set; }                           /// Путь к задаче 

    public string GetName();                                   /// Возвращает имя метода
    public Vector solve(bool log = false);                     /// Основная функция решения
}
