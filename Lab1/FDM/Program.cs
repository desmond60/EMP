try
{
    string json = File.ReadAllText(@"test\easy.json");          //: Простейший тест (1)
    //string json = File.ReadAllText(@"test\sec_kraev.json");     //: Тест на вторые краевые (2)
    //string json = File.ReadAllText(@"test\polynom_2.json");     //: Тест на полином второй степени (3)
    //string json = File.ReadAllText(@"test\polynom_3.json");     //: Тест на полином третьей степени (4)
    //string json = File.ReadAllText(@"test\polynom_4.json");     //: Тест на полином четвертой степени (5)
    //string json = File.ReadAllText(@"test\not_polynomh1.json"); //: Тест на не полиноминальной функции (h1 - шаг) (6)
    //string json = File.ReadAllText(@"test\not_polynomh2.json"); //: Тест на не полиноминальной функции (h1/2 - шаг) (6)
    //string json = File.ReadAllText(@"test\not_polynomh3.json"); //: Тест на не полиноминальной функции (h1/4 - шаг) (6)
    //string json = File.ReadAllText(@"test\uneven.json");        //: Тест на неравномерной сетке (1)

    Data data = JsonConvert.DeserializeObject<Data>(json)!;

    if (data is null) throw new FileNotFoundException("File uncorrected!");

    Solve task = new Solve(data, 1);
    task.solve();
}
catch (FileNotFoundException ex)
{
    WriteLine(ex.Message);
}
