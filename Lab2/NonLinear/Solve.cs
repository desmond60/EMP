public partial class Solve
{
    protected Data    data;                    /// Данные задачи
    protected Method  method;                  /// Метод решения (Итерация или Ньютон)
    protected Elem[]  elems;                   /// Конечные элементы
    protected SLAU    slau;                    /// Структура СЛАУ
    protected string  path  { get; set; }      /// Путь к папке с задачей

    public Solve(Data _data, Method _method, uint _funcNumber) {
        Function.NumberFunc = _funcNumber;
        method              = _method    ;
        data                = _data      ;
    }

    public partial void solve() {
        elems = Generate();              //? Генерация сетки
        memory();                        //? Выделение памяти под матрица и вектора

        IMethods task = method switch 
                    {
                        Method.Iteration => new Iteration(slau, data, elems),
                        Method.Newton    => new Newton(slau, data, elems),
                        _                => new Iteration(slau, data, elems)
                    };
        Directory.CreateDirectory($"{path}\\output");
        File.WriteAllText($"{path}\\output\\table_" + $"{task.GetName()}" + ".txt", task.solve());
    }

    private partial void memory() {
        slau.di   = new double[data.countEl + 1];
        slau.dl   = new double[data.countEl + 1];
        slau.du   = new double[data.countEl + 1];
        slau.q    = new double[data.countEl + 1];
        slau.q1   = new double[data.countEl + 1];
        slau.f    = new double[data.countEl + 1];
        // Записываем начальное приближение
        Array.Copy(data.init_q, slau.q, slau.q.Length);
    }
}
