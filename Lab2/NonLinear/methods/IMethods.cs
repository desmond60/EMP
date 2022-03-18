namespace NonLinear.methods;
public interface IMethods
{
    public string GetName();                                              /// Возвращает имя метода 
    public string solve();                                                /// Главная функция решения
    public void global();                                                 /// Генерация глобальной матрицы
    public (double[][], double[]) local(int _i);                          /// Генерация локальной матрицы
    public void add_to_global(double[][] _mat, double[] _vec,  int _i);   /// Занесение локальной в глоабльную 
    public void kraev();                                                  /// Учет краевых условий
    public void firstKraev(int _i);                                       /// Учет первого краевого
    public void secondKraev(int _i);                                      /// Учет второго краевого
    public void thirdKraev(int _i);                                       /// Учет третьего краевого
    public double[] newApprox();                                          /// Получение нового приближения (релаксация)
    public StringBuilder output();                                        /// Вывод решения
}