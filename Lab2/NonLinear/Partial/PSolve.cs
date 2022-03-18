public partial class Solve
{
    //* Генерация сетки
    private Elem[] Generate() {
        Elem[] elems = new Elem[data.countEl];

        // Подсчет начального шага
        double _h = data.k == 1 
                    ? (data.end - data.begin) / data.countEl
                    : (1 - data.k) * (data.end - data.begin) / (1 - Pow(data.k, data.countEl));

        double x = data.begin; // Текущая позиция
        for (int i = 0; i < data.countEl && x <= data.end; i++) {
            elems[i] = new Elem(x, x + _h, _h);
            x  += _h;
            _h *= data.k;
        }
        
        StringBuilder grid = new StringBuilder();
        grid.Append($"{"x1",20} {"x2",24} {"h",24}\n");
        grid.Append(String.Join("\n", elems));
        File.WriteAllText($"{path}\\grid.txt", grid.ToString());
        return elems;
    }

    //* Путь к папке, где наход. решение задачи
    public void SetPath(string _path) {
        path = _path;
    }

    public partial void solve();
    private partial void memory();
}