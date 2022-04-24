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
    Solve task = new Solve(grid);
    task.solve();

    // Solve task = args[3] switch
    // {
    //     "iteration" => new Solve(data, Method.Iteration, uint.Parse(args[5])),
    //     "newton"    => new Solve(data, Method.Newton,    uint.Parse(args[5])),
    //     _           => new Solve(data, Method.Iteration, uint.Parse(args[5]))
    // };
    // task.SetPath(Path.GetDirectoryName(args[1])!);
    // task.solve();
}
catch (FileNotFoundException ex) {
    WriteLine(ex.Message);
}
catch (ArgumentException ex) {
    ShowHelp();
    WriteLine(ex.Message);
}
