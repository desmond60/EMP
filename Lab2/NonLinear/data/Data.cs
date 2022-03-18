namespace NonLinear;
public class Data
{
    //* Данные задачи
    public uint     countEl { get; set; }   /// Количество КЭ/отрезков
    public double   begin   { get; set; }   /// Начало отрезка
    public double   end     { get; set; }   /// Конец отрезка
    public double   k       { get; set; }   /// Коэффициент разрядки
    public double   gamma   { get; set; }   /// Коэффициенты гамма
    public double   betta   { get; set; }   /// Коэффициенты бетта
    public double[] kraev   { get; set; }   /// Краевые условия
    public double[] init_q  { get; set; }   /// Начальное приближение
}