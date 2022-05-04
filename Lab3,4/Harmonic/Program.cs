using System.Security.Principal;
try {

    if (args.Length == 0) throw new ArgumentException("Not found arguments!");
    if (args[0] == "-help") {
        ShowHelp(); return;
    }

    string json = File.ReadAllText(args[1]);
    Data data = JsonConvert.DeserializeObject<Data>(json)!;
    if (data is null) throw new FileNotFoundException("File uncorrected!");

    // Определение функции
    Function.Init(data.N);

    // Генерация сетки
    Generate generator = new Generate(data, Path.GetDirectoryName(args[1])!);
    Grid grid = generator.generate();
    
    // Трансформация сетки (под синус и косинус) //: Для удобства
    grid = generator.transformation(grid);

    // Решение задачи
    Solve task = new Solve(grid, Path.GetDirectoryName(args[1])!);
    switch(args[2]) 
    {
        case "--los"     : task.solve(Method.LOS    ); break;
        case "--lu"      : task.solve(Method.LU     ); break;
        case "--bcgstab" : task.solve(Method.BCGSTAB); break;
        default          : task.solve(Method.LOS    ); break;
    };
}
catch (FileNotFoundException ex) {
    WriteLine(ex.Message);
}
catch (ArgumentException ex) {
    ShowHelp();
    WriteLine(ex.Message);
}
