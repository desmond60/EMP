namespace FDM;
public class Solve
{
    public int         CountX { get; set; }   /// Количество точек на оси X
    public int         CountY { get; set; }   /// Количество точек на оси Y
    public double[]    X      { get; set; }   /// Значения X-ов
    public double[]    Y      { get; set; }   /// Значения Y-ов
    public int         GX     { get; set; }   /// К-во точек на нижней границе области "Г"
    public int         GY     { get; set; }   /// К-во точек на правой границе области "Г" 
    private Matrix matrix;   /// 5-диагональная матрица 

    public Solve(Data data, uint Num) { 
        (CountX, CountY, X, Y, GX, GY) = data; 
        Function.NumberFunc = Num;
        Function.Init();
    }

    //* Решение задачи
    public void solve() {
        memory();                                        //? Выделение памяти
        completion();                                    //? Заполнение матрицы 
        var task = new Seidel(matrix, 10000, 1e-14);     //? Создание метода Гаусса-Зейделя
        task.solve(true);                                //? Решение СЛАУ
        writeTable();                                    //? Записб таблички с решением
    }

    //* Заполнение матрицы (область "Г") снизу->вверх
    private void completion() {
        int id = 0;                //: индекс узла
        double hx1, hx2, hy1, hy2; //: h-ки приращения аргументов
        
        // Нижняя линия области "Г"
        matrix.pr[0] = Absolut(X[0], Y[0]);       // Левый нижний угол
        hy1 = Abs(Y[0] - Y[1]);
        for (int i = 1; i < GX - 1; i++) {
            matrix.pr[i] = Absolut(X[i], Y[0], 1);
            if (IsSecondKraev(1)) {
                matrix.di[i]  = -lambda / hy1;
                matrix.du2[i] =  lambda / hy1; 
            }
        }
        id = GX - 1;
        matrix.pr[id] = Absolut(X[GX - 1], Y[0]);  // Правый нижний угол
        id++;

        // До линии между шапкой и ножкой
        for (int i = 1; i < GY - 1; i++, id++) {
            hy1 = Abs(Y[i] - Y[i - 1]);
            hy2 = Abs(Y[i + 1] - Y[i]);
            matrix.pr[id] = Absolut(X[0], Y[i], 2);
            if (IsSecondKraev(2)) {
                hx1 = Abs(X[0] - X[1]);
                matrix.di [id] = -lambda / hx1;
                matrix.du1[id] =  lambda / hx1;
            }
            id++;

            for (int j = 1; j < GX - 1; j++, id++) {
                hx1 = Abs(X[j] - X[j - 1]);
                hx2 = Abs(X[j + 1] - X[j]);
                matrix.pr [id]                 = Func(X[j], Y[i]);
                matrix.dl1[id - 1]             = -2*lambda / (hx1 * (hx2 + hx1)); 
                matrix.dl2[id - matrix.shift1] = -2*lambda / (hy1 * (hy2 + hy1));
                matrix.du1[id]                 = -2*lambda / (hx2 * (hx2 + hx1));
                matrix.du2[id]                 = -2*lambda / (hy2 * (hy2 + hy1));
                matrix.di [id]                 = lambda * (2/(hx1*hx2) + 2/(hy1*hy2)) + gamma; 
            }

            matrix.pr[id] = Absolut(X[GX - 1], Y[i], 3);
            if (IsSecondKraev(3)) {
                hx1 = Abs(X[GX - 1] - X[GX - 2]);
                matrix.di[id] = lambda / hx1;
                matrix.dl1[id - 1] = -lambda / hx1;
            }
        }
    
        matrix.dshift = id;

        // Между шляпкой и ножкой 
        hy1 = Abs(Y[GY - 1] - Y[GY - 2]);
        hy2 = Abs(Y[GY] - Y[GY - 1]);
        matrix.pr[id] = Absolut(X[0], Y[GY - 1], 2);
        if (IsSecondKraev(2)) {
            hx1 = Abs(X[0] - X[1]);
            matrix.di[id] = -lambda / hx1;
            matrix.du1[id] = lambda / hx1;
        }
        id++;

        for (int i = 1; i < GX; i++, id++) {
            hx1 = Abs(X[i] - X[i - 1]);
            hx2 = Abs(X[i + 1] - X[i]);
            matrix.pr[id]                  = Func(X[i], Y[GY - 1]);
            matrix.dl1[id - 1]             = -2*lambda / (hx1 * (hx2 + hx1)); 
            matrix.dl2[id - matrix.shift1] = -2*lambda / (hy1 * (hy2 + hy1));
            matrix.du1[id]                 = -2*lambda / (hx2 * (hx2 + hx1));
            matrix.du2[id]                 = -2*lambda / (hy2 * (hy2 + hy1));
            matrix.di [id]                 = lambda * (2/(hx1*hx2) + 2/(hy1*hy2)) + gamma; 
        }

        for (int i = GX; i < CountX - 1; i++, id++) {
            matrix.pr[id] = Absolut(X[i], Y[GY - 1], 4);
            if (IsSecondKraev(4)) {
                matrix.di[id] = -lambda / hy2;
                matrix.du2[id] = lambda / hy2;
            }
        }
        matrix.pr[id] = Absolut(X[CountX - 1], Y[GY - 1]);
        id++;

        // Шляпка 
        for (int i = GY; i < CountY - 1; i++) {
            hy1 = Abs(Y[i] - Y[i - 1]);
            hy2 = Abs(Y[i + 1] - Y[i]);
            matrix.pr[id] = Absolut(X[0], Y[i], 2);
            if (IsSecondKraev(2)) {
                hx1 = Abs(X[0] - X[1]);
                matrix.di[id] = -lambda / hx1;
                matrix.du1[id] = lambda / hx1;
            }    
            id++;
        
            for (int j = 1; j < CountX - 1; j++, id++) {
                hx1 = Abs(X[j] - X[j - 1]);
                hx2 = Abs(X[j + 1] - X[j]);
                matrix.pr[id]                  = Func(X[j], Y[i]);
                matrix.dl1[id - 1]             = -2*lambda / (hx1 * (hx2 + hx1)); 
                matrix.dl2[id - matrix.shift1] = -2*lambda / (hy1 * (hy2 + hy1));
                matrix.du1[id]                 = -2*lambda / (hx2 * (hx2 + hx1));
                matrix.du2[id]                 = -2*lambda / (hy2 * (hy2 + hy1));
                matrix.di [id]                 = lambda * (2/(hx1*hx2) + 2/(hy1*hy2)) + gamma; 
            }

            matrix.pr[id] = Absolut(X[CountX - 1], Y[i], 5);
            if (IsSecondKraev(5)) {
                hx1 = Abs(X[CountX - 1] - X[CountX - 2]);
                matrix.di[id] = lambda / hx1;
                matrix.dl1[id - 1] = -lambda / hx1;
            }
            id++;
        }

        // Верхушка шляпки
        matrix.pr[id] = Absolut(X[0], Y[CountY - 1]);
        id++;
        hy1 = Abs(Y[CountY - 1] - Y[CountY - 2]);
        for (int i = 1; i < CountX - 1; i++, id++) {
            matrix.pr[id] = Absolut(X[i], Y[CountY - 1], 6);
            if (IsSecondKraev(6)) {
                matrix.di[id] = lambda / hy1;
                matrix.dl2[id - matrix.shift1] = -lambda / hy1;
            }
        }
        matrix.pr[id] = Absolut(X[CountX - 1], Y[CountY - 1]);
    }

    //* Выделяем память под матрицу
    private void memory() {
        matrix.N = GX * (CountY - GY) + CountX * GY;  // Размерность матрицы
        matrix.shift1    = GX;
        matrix.shift2    = CountX;
        matrix.di        = new double[matrix.N];
        matrix.du1       = new double[matrix.N];
        matrix.du2       = new double[matrix.N];
        matrix.dl1       = new double[matrix.N];
        matrix.dl2       = new double[matrix.N];
        matrix.pr        = new double[matrix.N];
        matrix.x         = new double[matrix.N];
        matrix.absolut_x = new double[matrix.N];
        Array.Fill(matrix.di, 1); // Заполнение диагонали единицами
    }

    //* Заполнение и запись таблички с решением
    private void writeTable() {
        StringBuilder table = new StringBuilder();
        string margin = String.Join("", Enumerable.Repeat("-", 16));

        table.Append(String.Join("", Enumerable.Repeat("-", 86)) + "\n");
        table.Append($"|X{" ",-14} | Y{" ",-13} | U{" ",-12}  | U`{" ",-12} | |U`- U| {" ",-7}|\n");
        table.Append($"|" + margin + "|" + margin + "|" + margin + "|" + margin + "|" + margin + "|\n");

        int id      = 0;
        for (int i = 0; i < CountY - GY; i++) {
            for (int j = 0; j < GX; j++, id++) {
                double absolut = Absolut(X[j], Y[i]);
                table.Append($"|{String.Format("{0,16}", X[j])}" +
                                $"|{String.Format("{0,16}", Y[i])}" +
                                $"|{String.Format("{0,16}", matrix.x[id].ToString("E6"))}" + 
                                $"|{String.Format("{0,16}", absolut.ToString("E6"))}" + 
                                $"|{String.Format("{0,16}", Abs(absolut - matrix.x[id]).ToString("E6"))}|\n");
            }
        }

        for (int i = GY - 1; i < CountY; i++) {
            for (int j = 0; j < CountX; j++, id++) {
                double absolut = Absolut(X[j], Y[i]);
                table.Append($"|{String.Format("{0,16}", X[j])}" +
                                $"|{String.Format("{0,16}", Y[i])}" +
                                $"|{String.Format("{0,16}", matrix.x[id].ToString("E6"))}" + 
                                $"|{String.Format("{0,16}", absolut.ToString("E6"))}" + 
                                $"|{String.Format("{0,16}", Abs(absolut - matrix.x[id]).ToString("E6"))}|\n");
            }
        }
        table.Append(String.Join("", Enumerable.Repeat("-", 86)) + "\n");
        File.WriteAllText("test/tables/table.txt", table.ToString());
    }

}