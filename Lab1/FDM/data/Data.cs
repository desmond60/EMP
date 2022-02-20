namespace FDM;
public class Data
{
    //* Данные задачи
    public int         CountX { get; set; }   /// Количество точек на оси X
    public int         CountY { get; set; }   /// Количество точек на оси Y
    public double[]    X      { get; set; }   /// Значения X-ов
    public double[]    Y      { get; set; }   /// Значения Y-ов
    public int         GX     { get; set; }   /// К-во точек на нижней границе области "Г"
    public int         GY     { get; set; }   /// К-во точек на правой границе области "Г" 

    public void Deconstruct(out int      countX, 
                            out int      countY, 
                            out double[] x,
                            out double[] y, 
                            out int      gx,
                            out int      gy) 
    {
        countX = CountX;
        countY = CountY;
        x = X;
        y = Y;
        gx = GX;
        gy = GY;
    }
}