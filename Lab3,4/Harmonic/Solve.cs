using System.Runtime.Serialization;
using System.Diagnostics;
namespace Harmonic;

public class Solve
{
    protected SLAU slau;                                      /// Структура СЛАУ
    protected Grid grid;                                      /// Сетка (узлы, КЭ, краевые)

    private string Path;                                      /// Путь к задаче

    private Stopwatch time;                                   /// Для замера времени

    //: Табличка с окончательным решением
    StringBuilder table_sol;                                  /// Табличка с окончательным решением
    StringBuilder table_rate;                                 /// Табличка с относительной погрешностью
    //: ----------------------------------------------------------

    //* Конструктор
    public Solve(Grid grid, string path) {
        this.grid = grid;
        this.Path = path;

        time = new Stopwatch();
    }

    //* Основной метод решения
    public void solve(Method method) {
        portrait();           //? Составление портрета матрицы
        global();             //? Составление глобальной матрицы

        // Новый путь к решению
        Path = method switch
        {
            Method.LOS     => Path + "/output_LOS",
            Method.LU      => Path + "/output_LU",
            Method.BCGSTAB => Path + "/output_BCGSTAB"
        };
        Directory.CreateDirectory(Path);  //? Создание директории с решением

        IMethods issue = method switch    //? Определение метода
        {
            Method.LOS     => new LOS(slau, 1e-16, 10000),
            Method.LU      => new LU(slau),
            Method.BCGSTAB => new BCG(slau, 1e-16, 30000),
            _ => new LOS(slau, 1e-16, 10000)
        };
        issue.Path = this.Path;            //? Установление пути
        time.Start();                      //? Запуск таймера
        slau.q = issue.solve(true);        //? Решение (true - записать итерации)
        time.Stop();                       //? Остановка таймера
        output();                          //? Запись в файл решение
    }

    //* Составление портрета матрицы (ig, jg, выделение памяти)
    private void portrait() {
        Portrait port = new Portrait(grid.Count_Node);

        // Генерируем массивы ig и jg и размерность
        slau.N_el = port.GenPortrait(ref slau.ig, ref slau.jg, grid.Elems);
        slau.N    = grid.Count_Node;

        // Выделяем память
        slau.gl = new Vector(slau.N_el);
        slau.gu = new Vector(slau.N_el);
        slau.di = new Vector(slau.N);
        slau.f  = new Vector(slau.N);
        slau.q  = new Vector(slau.N);
    }

    //* Составление глобальной матрицы
    private void global() {

        // Для каждого конечного элемента
        for (int index_fin_el = 0; index_fin_el < grid.Count_Elem; index_fin_el++) {
            // Составление локальной матрицы
            (double[][] local_matrix, Vector local_f) = local(index_fin_el);

            // Добавление в глобальную
            EntryGlobalMatrix(local_matrix, local_f, index_fin_el);
        }

        // Учет только первых краевых 
        for (int index_kraev = 0; index_kraev < grid.Count_Kraev; index_kraev++)
            FirstKraev(index_kraev);   
    }

    //* Составление локальной матрицы и локального вектора
    private (double[][], Vector) local(int index) {
        
        // Локальные функции
        int mu(int i) => i % 2;
        int v (int i) => (i/4) % 2;
        int up(int i) => (i/2) % 2;
        
        double[][] mat = new double[16][];
        mat = mat.Select(n => new double[16]).ToArray();
        Vector vec = new Vector(16);

        // Подсчет lambda, sigma, hi
        double lambda = 0, sigma = 0, hi = 0;
        Vector f_sin = new Vector(8);
        Vector f_cos = new Vector(8);
        for (int i = 0; i < 8; i++) {
            double x = grid.Nodes[grid.Elems[index].Node[2*i]].x;
            double y = grid.Nodes[grid.Elems[index].Node[2*i]].y;
            double z = grid.Nodes[grid.Elems[index].Node[2*i]].z;
            var point = new Vector(new double[]{x, y, z});
            lambda += Lambda(point);
            sigma  += Sigma (point);
            hi     += Hi    (point);

            // Подсчет f
            f_sin[i] = Fs(point);
            f_cos[i] = Fc(point);
        }
        lambda /= 8.0; sigma /= 8.0; hi /= 8.0;

        // Подсчет hx, hy, hz
        double hx = Abs(grid.Nodes[grid.Elems[index].Node[2]].x - grid.Nodes[grid.Elems[index].Node[0]].x);
        double hy = Abs(grid.Nodes[grid.Elems[index].Node[4]].y - grid.Nodes[grid.Elems[index].Node[0]].y);
        double hz = Abs(grid.Nodes[grid.Elems[index].Node[8]].z - grid.Nodes[grid.Elems[index].Node[0]].z);

        // Построение матрицы жесткости и массы        
        double[,] G1 = new double[2,2] {{1, -1}, {-1, 1} };
        double[,] M1 = new double[2,2] {{1/3.0, 1/6.0}, {1/6.0, 1/3.0} };
        
        double[][] G = new double[8][];
        G = G.Select(n => new double[8]).ToArray();
        double[][] M = new double[8][];
        M = M.Select(n => new double[8]).ToArray();

        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++) {
                G[i][j]  = (hy*hz/hx) * G1[mu(i), mu(j)] * M1[v(i), v(j)] * M1[up(i), up(j)];
                G[i][j] += (hx*hz/hy) * M1[mu(i), mu(j)] * G1[v(i), v(j)] * M1[up(i), up(j)];
                G[i][j] += (hy*hx/hz) * M1[mu(i), mu(j)] * M1[v(i), v(j)] * G1[up(i), up(j)];           
                G[i][j] *= lambda;
                M[i][j]  = hx*hy*hz * M1[mu(i), mu(j)] * M1[v(i), v(j)] * M1[up(i), up(j)];
            }

        // Построение локальной матрицы
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++) {
                mat[2*i + 1][2*j + 1] = mat[2*i][2*j] = G[i][j] - omega*omega*hi*M[i][j];
                mat[  2*i  ][2*j + 1] = -omega*sigma*M[i][j];
                mat[2*i + 1][  2*j  ] = omega*sigma*M[i][j];
            }

        // Подсчет правой части
        Vector b_sin = new Vector(8);
        Vector b_cos = new Vector(8);
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++) {
                b_sin[i] += M[i][j]*f_sin[j];
                b_cos[i] += M[i][j]*f_cos[j];
            }
        
        // Построение правой части
        for (int i = 0; i < 8; i++)
            (vec[2*i], vec[2*i + 1]) = (b_sin[i], b_cos[i]);

        return (mat, vec);
    }

    //* Поиск столбца
    private int find(int f, int s) {
        int col = 0;
        for (int i = slau.ig[f]; i < slau.ig[f + 1]; i++) {
            if (slau.jg[i] == s) {
                col = i;
                break;
            }
        }
        return col;
    }

    //* Занесение матрицы и вектора в глобальную
    private void EntryGlobalMatrix(double[][] mat, Vector vec, int index) {
        for (int i = 0; i < 16; i++) {
            int row = grid.Elems[index].Node[i];
            for (int j = 0; j < i; j++) {
                if (row > grid.Elems[index].Node[j]) {
                    int col = find(row, grid.Elems[index].Node[j]);
                    slau.gl[col] += mat[i][j];
                    slau.gu[col] += mat[j][i];
                }
                else {
                    int col = find(grid.Elems[index].Node[j], row);
                    slau.gu[col] += mat[i][j];
                    slau.gl[col] += mat[j][i];
                }
            }
            slau.di[row] += mat[i][i];
            slau.f[row]  += vec[i];
        }
    }

    //* Учет первых краевых 
    private void FirstKraev(int index) {
        (int row, double value) = (grid.Kraevs[index].Node[0], grid.Kraevs[index].Value);

        // Стави значения
        slau.di[row] = 1;
        slau.f [row] = value;
        
        // Зануляем в треугольнике (столбца)
        for (int i = slau.ig[row]; i < slau.ig[row + 1]; i++) {
            slau.f[slau.jg[i]] -= slau.gu[i]*value;
            slau.gl[i] = 0;
            slau.gu[i] = 0;
        }

        // Зануляем в треугольнике (строки)
        for (int i = row + 1; i < slau.N; i++) {
            for (int j = slau.ig[i]; j < slau.ig[i + 1]; j++) {
                if (slau.jg[j] == row) {
                    slau.f[i] -= slau.gl[j]*value;
                    slau.gl[j] = 0;
                    slau.gu[j] = 0;
                }
            }
        }
    }   

    //* Запись решения
    private void output() {
        //: Табличка с окончательным решением
        string margin = String.Join("", Enumerable.Repeat("-", 16));
        
        table_sol = new StringBuilder($"|U{" ",-14} | U`{" ",-12} | |U`- U| {" ",-7}|\n");
        table_sol.Append($"|" + margin + "|" + margin + "|" + margin + "|\n");

        
        for (int i = 0; i < slau.N; i++) {
            Vector node = new Vector(new double[] { grid.Nodes[i].x, grid.Nodes[i].y, grid.Nodes[i].z });
            if (i % 2 == 0) {
                table_sol.Append($"|{String.Format("{0,16}", slau.q[i].ToString("E6"))}"  + 
                                 $"|{String.Format("{0,16}", Us(node).ToString("E6"))}"   + 
                                 $"|{String.Format("{0,16}", Abs(Us(node) - slau.q[i]).ToString("E6"))}|\n");
            } else {
                table_sol.Append($"|{String.Format("{0,16}", slau.q[i].ToString("E6"))}"  + 
                                 $"|{String.Format("{0,16}", Uc(node).ToString("E6"))}"   + 
                                 $"|{String.Format("{0,16}", Abs(Uc(node) - slau.q[i]).ToString("E6"))}|\n");                
            }
        }

        table_sol.Append(String.Join("", Enumerable.Repeat("-", 52)) + "\n");

        File.WriteAllText(Path + "/table_sol.txt", table_sol.ToString());
        //: ---------------------------------------------------------

        //: Табличка с погрешностью 
        table_rate = new StringBuilder();

        double diffetence_sin, difference_cos, rate_u_sin, rate_u_cos;
        difference_cos = diffetence_sin = rate_u_cos = rate_u_sin = 0;

        for (int i = 0; i < grid.Count_Node / 2; i++) {
            Vector node = new Vector(new double[] { grid.Nodes[2*i].x, grid.Nodes[2*i].y, grid.Nodes[2*i].z });
            double us = Us(node);
            double uc = Uc(node);
            diffetence_sin += Pow((us - slau.q[2*i]), 2);
            difference_cos += Pow((uc - slau.q[2*i + 1]), 2);
            rate_u_sin     += Pow(us, 2);
            rate_u_cos     += Pow(uc, 2);
        }

        table_rate.Append($"Общее:   {Sqrt((diffetence_sin + difference_cos) / (rate_u_sin + rate_u_cos)).ToString("E6")}\n");
        table_rate.Append($"Синус:   {Sqrt(diffetence_sin / rate_u_sin).ToString("E6")}\n");
        table_rate.Append($"Косинус: {Sqrt(difference_cos / rate_u_cos).ToString("E6")}\n");
        
        TimeSpan ts = time.Elapsed;
        string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
            ts.Hours, ts.Minutes, ts.Seconds,
            ts.Milliseconds / 10);
        table_rate.Append($"Время: {elapsedTime}\n");

        File.WriteAllText(Path + "/table_rate.txt", table_rate.ToString());
        //: ----------------------------------------------------------
    }
}
