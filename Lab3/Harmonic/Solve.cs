namespace Harmonic;

public class Solve
{
    SLAU slau;                                      /// Структура СЛАУ
    Grid grid;                                      /// Сетка (узлы, КЭ, краевые)

    //* Конструктор
    public Solve(Grid grid) {
        this.grid = grid;
    }

    //* Основной метод решения
    public void solve() {
        portrait();           //? Составление портрета матрицы
    }

    //* Составление портрета матрицы (ig, jg, выделение памяти)
    private void portrait() {
        Portrait port = new Portrait(grid.Count_Node);

        // Генерируем массивы ig и jg и размерность
        slau.N = port.GenPortrait(ref slau.ig, ref slau.jg, grid.Elems);
    }


}
