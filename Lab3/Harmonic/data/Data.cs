namespace Harmonic;
public class Data
{
    //* Данные для генерации сетки
    public double[] start { get; set; }       /// Начальная точка
    public double[] end   { get; set; }       /// Конечная точка
    public double   hx    { get; set; }       /// Шаг по Оси X
    public double   hy    { get; set; }       /// Шаг по Оси Y
    public double   hz    { get; set; }       /// Шаг по Оси Z 
    public double   kx    { get; set; }       /// Коэффициент деления по Оси X  
    public double   ky    { get; set; }       /// Коэффициент деления по Оси Y  
    public double   kz    { get; set; }       /// Коэффициент деления по Оси Z 

    //* Деконструктор
    public void Deconstruct(out Vector   Start, 
                            out Vector   End,
                            out double   Hx, 
                            out double   Hy, 
                            out double   Hz, 
                            out double   Kx,
                            out double   Ky,
                            out double   Kz) 
    {
        Start = new Vector(start);
        End   = new Vector(end);
        Hx    = hx;
        Hy    = hy;
        Hz    = hz;
        Kx    = kx;
        Ky    = ky;
        Kz    = kz;
    }
}