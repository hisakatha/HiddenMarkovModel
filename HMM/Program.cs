using System;
using System.Linq;
using System.Collections.Generic;
//using System.Text;

//using System.Array;
// 以下で使われるSystem.Array.IndexOf(array, value)メソッドは、
// 配列array中に値valueが最初に現れるインデックスを返すメソッドであり、名前空間System.Arrayに属する。

namespace HMM
{
    /// <summary>
    /// Mainメソッドを持ちファイル読み込み等を行うクラス。
    /// </summary>
    class Program
    {
        // エラー時に表示する利用方法の文字列。
        static readonly string usage = "Usage: command {parameter file} {FASTA file(s)}...";

        /// <summary>
        /// 1つ目のファイルパスからHMMのパラメータを読み取り、HMMを作成し、2つ目のfasta形式のファイルパスから観測文字列を読み取る。
        /// このHMMに観測文字列を与えることで、隠れ状態の推定などを行う。
        /// </summary>
        /// <param name="args">コマンドライン引数</param>
        static void Main(string[] args)
        {
            if (args.Length < 2)
            {
                Console.WriteLine("Wrong usage.");
                Console.WriteLine(usage);
                return;
            }
            // HMMクラスの変数を宣言。
            HMM model;
            int index = 0;
            int blockLength = 1;
            
            // try-catch文
            // 例外発生時にcatch文が実行される。
            try
            {
                // outキーワードは参照渡しを指示するキーワードであるが、用途が限定されている。
                // outキーワードが付いた変数は呼び出し元で初期化する必要がなく、必ずメソッド側で値が設定されることが約束される。
                ParseParameter(out model, args[index]);
                index++;
            }
            catch
            {
                // モデル作成に失敗した場合の処理
                return;
            }

            string[] sequenceArray;
            string[] fastas = new string[args.Length - index];
            // 配列fastasに、コマンドライン引数から該当するFASTAファイルパスをコピーする。
            System.Array.Copy(args, index, fastas, 0, args.Length - index);
            // try-catch文
            // 例外発生時にcatch文が実行される。
            try
            {
                sequenceArray = ParseFastaFiles(fastas);
            }
            catch
            {
                return;
            }
            
            // 各塩基配列についてブロック長を大きくしながら複数回時間計測を行う。
            // 1条件での繰り返し数
            const int repeatNumber = 10;
            // 通常のViterbiアルゴリズムでの実行時間を格納する配列。
            long[] normalTimes = new long[sequenceArray.Length];
            // 圧縮を用いたViterbiアルゴリズムでの実行時間を格納するQueue。
            Queue<long>[] times = new Queue<long>[sequenceArray.Length];
            for (int i = 0; i < times.Length; i++)
            {
                times[i] = new Queue<long>();
            }
            // .NETに用意されているStopwatchクラスのインスタンスを生成する。
            System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
            for (int i = 0; i < sequenceArray.Length; i++)
            {
                timer.Reset();
                for (int t = 0; t < repeatNumber; t++)
                {
                    timer.Start();
                    model.Viterbi(sequenceArray[i]);
                    timer.Stop();
                }
                normalTimes[i] = timer.ElapsedMilliseconds;
                
                int maxBlockLength = (int)System.Math.Ceiling(System.Math.Log(sequenceArray[i].Length, model.AlphabetSize));
                for (blockLength = 2; blockLength <= maxBlockLength; blockLength++)
                {
                    timer.Reset();
                    for (int t = 0; t < repeatNumber; t++)
                    {
                        timer.Start();
                        model.ViterbiWithCompression(sequenceArray[i], blockLength);
                        timer.Stop();
                    }
                    times[i].Enqueue(timer.ElapsedMilliseconds);
                }
            }

            while (times.Any((Queue<long> element) => (element.Count > 0)))
            {
                for (int i = 0; i < sequenceArray.Length; i++)
                {
                    if (times[i].Count > 0)
                    {
                        Console.Write("{0:G5}", (double)normalTimes[i] / (double)times[i].Dequeue());
                    }
                    else
                    {
                        Console.Write("N/A");
                    }
                    Console.Write("\t");
                }
                Console.WriteLine("");
            }

#if DEBUG
            Console.ReadKey(true);
#endif
            return;
        }

        /// <summary>
        /// 配列の要素を半角スペース区切りで表示するジェネリックメソッド。
        /// </summary>
        public static void PrintArray<type>(type[] array)
        {
            // バッファーサイズの初期値を1要素当たり4バイトとする。
            System.Text.StringBuilder buf = new System.Text.StringBuilder(array.Length * 4);
            for (int i = 0; i < array.Length - 1; i++)
            {
                buf.Append(array[i].ToString()).Append(" ");
            }
            buf.Append(array[array.Length - 1].ToString());
            Console.WriteLine(buf.ToString());
            return;
        }

        /// <summary>
        /// 二次元配列の要素を空白区切りで1行ずつ表示するジェネリックメソッド
        /// </summary>
        /// <typeparam name="type">二次元配列の型</typeparam>
        /// <param name="matrix">要素を表示する二次元配列</param>
        public static void PrintMatrix<type>(type[][] matrix)
        {
            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j < matrix[i].Length; j++)
                {
                    Console.Write("{0:F10} ", matrix[i][j]);
                }
                Console.WriteLine("");
            }
            return;
        }

        /// <summary>
        /// パスで指定したファイルからHMMのパラメータを読み取り、そのパラメータでHMMの実体（インスタンス）を生成する。
        /// </summary>
        /// <param name="model">HMMクラスのインスタンスが格納される変数</param>
        /// <param name="path">パラメータのファイルパス</param>
        static void ParseParameter(out HMM model, string path)
        {
            // ファイル読み出しのストリームを用意。
            System.IO.StreamReader stream = null;
            // try-catch-finally文
            // 例外発生時にcatch文が実行され、例外の有無に関わらずfinally文が実行される。
            try
            {
                // 読み込み専用で"path"で指定されたファイルを開く。
                stream = new System.IO.StreamReader(path);
                // bufferArrayにアルファベットサイズ、アルファベット、状態数を示す文字列を格納する。
                string[] bufferArray = new string[3];
                for (int i = 0; i < 3; i++)
                {
                    while (true)
                    {
                        bufferArray[i] = stream.ReadLine().Trim();
                        int index = bufferArray[i].IndexOf('%');
                        if (index == -1)
                        {
                            break;
                        }
                        else if (index > 0)
                        {
                            bufferArray[i] = bufferArray[i].Remove(index).Trim();
                            break;
                        }
                    }
                }
                // 読み取った文字列を数値等に変換し、HMMクラスのインスタンスを生成する。
                model = new HMM(int.Parse(bufferArray[0]), bufferArray[1].Replace(" ", "").ToCharArray(), int.Parse(bufferArray[2]));
                // 遷移確率を読み取り、モデルに設定する。
                string buffer;
                for (int i = 0; i < model.StateNumber; i++)
                {
                    while (true)
                    {
                        buffer = stream.ReadLine().Trim();
                        int index = buffer.IndexOf('%');
                        if (index == -1)
                        {
                            break;
                        }
                        else if (index > 0)
                        {
                            buffer = buffer.Remove(index).Trim();
                            break;
                        }
                    }
                    string[] strArray = buffer.Split(' ');
                    if (strArray.Length == model.StateNumber)
                    {
                        for (int j = 0; j < model.StateNumber; j++)
                        {
                            model.SetTransitionProbability(i, j, double.Parse(strArray[j]));
                        }
                    }
                    else
                    {
                        Console.WriteLine("A dimension of transition probability matrix is wrong.");
                    }
                }
                // 出力確率を読み取り、モデルに設定する
                for (int i = 0; i < model.StateNumber - 1; i++)
                {
                    while (true)
                    {
                        buffer = stream.ReadLine().Trim();
                        int index = buffer.IndexOf('%');
                        if (index == -1)
                        {
                            break;
                        }
                        else if (index > 0)
                        {
                            buffer = buffer.Remove(index).Trim();
                            break;
                        }
                    }
                    string[] strArray = buffer.Split(' ');
                    if (strArray.Length == model.AlphabetSize)
                    {
                        for (int j = 0; j < model.AlphabetSize; j++)
                        {
                            model.SetEmissionProbability(j, i, double.Parse(strArray[j]));
                        }
                    }
                    else
                    {
                        Console.WriteLine("A dimension of emission probability matrix is wrong.");
                    }
                }
            }
            catch (Exception e)
            {
                // ファイル読み込みに失敗した場合の処理
                Console.WriteLine("Parameter file : {0}", e.Message);
                model = null;
                throw;
            }
            finally
            {
                // ファイル読み込みの成功に関わらず実行する処理
                if (stream != null)
                {
                    stream.Close();
                }
            }
            return;
        }

        /// <summary>
        /// FASTAファイルパスの格納された引数を受け取り、そこから配列を読み取り、string型の配列にして返す。
        /// </summary>
        /// <param name="args">ファイルパスの配列</param>
        /// <remarks>Ver.1.1 : 空行を読み飛ばすよう変更</remarks>
        static string[] ParseFastaFiles(string[] args)
        {
            // string型のリストを宣言。C++でのVectorと同じ。
            System.Collections.Generic.List<string> sequenceList = new System.Collections.Generic.List<string>();

            // ファイル読み出しのストリームを用意。
            System.IO.StreamReader stream = null;
            // 観測文字列をファイルから読み込む。
            // string型より文字列操作において高速なStringBuilderクラスを利用する。
            const int defaultCapacity = 100000;
            System.Text.StringBuilder seq = new System.Text.StringBuilder(defaultCapacity);
            // try-catch-finally文
            // 例外発生時にcatch文が実行され、例外の有無に関わらずfinally文が実行される。
            try
            {
                for (int i = 0; i < args.Length; i++)
                {
                    stream = new System.IO.StreamReader(args[i]);
                    // FASTA形式ファイルの最初の1行を読み飛ばす。
                    stream.ReadLine();
                    seq.Clear();
                    string buffer;
                    while (!stream.EndOfStream)
                    {
                        buffer = stream.ReadLine();
                        if (buffer.Length == 0)
                        {
                            // 空行を読み飛ばす。
                            continue;
                        }
                        if (buffer[0] == '>')
                        {
                            // '>'を行頭に見つけたら読み取った配列をsequenceListに追加し、新しく読み込みをseqに貯める。
                            sequenceList.Add(seq.ToString());
                            seq.Clear();
                        }
                        else
                        {
                            seq.Append(buffer);
                        }
                    }
                    stream.Close();
                    sequenceList.Add(seq.ToString());
                }
            }
            catch (Exception e)
            {
                // 例外処理
                Console.WriteLine("FASTA file : {0}", e.Message);
                throw;
            }
            finally
            {
                // 必ず行われる処理
                if (stream != null)
                {
                    stream.Close();
                }
            }
            // リストを通常の配列に変換して返す。
            return sequenceList.ToArray();
        }
    }

    /// <summary>
    /// HMMのパラメータを保有し、各種メソッドを提供するクラス
    /// </summary>
    class HMM
    {
        // readonlyキーワードによりコンストラクタ以降に変更できないメンバ変数を用意する。
        // public等を付加しない限り、メンバ変数はprivateである。
        readonly int alphabetSize;
        readonly char[] alphabet;
        readonly int stateNumber;
        double[][] transitionProbability;
        double[][] emissionProbability;

        /// <summary>
        /// コンストラクタ
        /// </summary>
        /// <param name="alphabetSize">アルファベットサイズ</param>
        /// <param name="alphabet">アルファベットのchar型配列</param>
        /// <param name="stateNumber">状態数</param>
        public HMM(int alphabetSize, char[] alphabet, int stateNumber)
        {
            this.alphabetSize = alphabetSize;
            this.stateNumber = stateNumber;
            this.alphabet = alphabet;
            // 確率行列の次元をここで決める。
            transitionProbability = new double[stateNumber][];
            for (int i = 0; i < stateNumber; i++)
            {
                transitionProbability[i] = new double[stateNumber];
            }
            emissionProbability = new double[alphabetSize][];
            for (int i = 0; i < alphabetSize; i++)
            {
                emissionProbability[i] = new double[stateNumber - 1];
            }
        }

        /// <summary>
        /// alphabetSizeに対するgetアクセサーを宣言し、外部からプロパティとして参照できるようにする。
        /// </summary>
        /// <remarks>
        /// C++でのgetter。
        /// 
        /// int AlphabetSize(void)
        /// {
        ///     return alphabetSize;
        /// }
        /// と同じだが、こちらの方が最適化されやすい。
        /// </remarks>
        public int AlphabetSize
        {
            get
            {
                return alphabetSize;
            }
        }

        /// <summary>
        /// stateNumberに対するgetアクセサーを宣言し、外部からプロパティとして参照できるようにする。
        /// </summary>
        public int StateNumber
        {
            get
            {
                return stateNumber;
            }
        }

        /// <summary>
        /// transitionProbabilityに対するgetアクセサーを宣言し、
        /// 外部からプロパティとして参照できるようにする。
        /// </summary>
        public double[][] TransitionProbability
        {
            get
            {
                return transitionProbability;
            }
        }

        /// <summary>
        /// emissionProbabilityに対するgetアクセサーを宣言し、外部からプロパティとして参照できるようにする。
        /// </summary>
        public double[][] EmissionProbability
        {
            get
            {
                return emissionProbability;
            }
        }

        /// <summary>
        /// 状態iから状態jへ遷移する確率をpにセットする。
        /// </summary>
        /// <param name="i">遷移前の状態</param>
        /// <param name="j">遷移後の状態</param>
        /// <param name="p">遷移確率</param>
        public void SetTransitionProbability(int i, int j, double p)
        {
            transitionProbability[i][j] = p;
            return;
        }

        /// <summary>
        /// x番目のアルファベットが状態kから出力される確率をpにセットする。
        /// </summary>
        /// <param name="x">出力されるアルファベットの0から始まるインデックス</param>
        /// <param name="k">出力時の状態</param>
        /// <param name="p">出力確率</param>
        public void SetEmissionProbability(int x, int k, double p)
        {
            emissionProbability[x][k] = p;
            return;
        }

        int SubstringToIndex(string sequence, int startIndex, int length)
        {
            int index = 0;
            for (int i = startIndex;  i < sequence.Length && i < startIndex + length; i++)
            {
                index = index * alphabetSize + System.Array.IndexOf<char>(alphabet, sequence[i]);
            }
            for (int i = sequence.Length; i < startIndex + length; i++)
            {
                index = index * alphabetSize;
            }
            return index;
        }

        static int PowerInt(int x, int y)
        {
            int result = 1;
            for (int i = 0; i < y; i++)
            {
                result *= x;
            }
            return result;
        }

        public int[] ViterbiWithCompression(string emission, int blockLength, bool debugmode = false)
        {
            // int patternNumber = PowerInt(alphabetSize, blockLength);
            // ブロックの作成
            double[][][] blockLogViterbi = new double[alphabetSize][][];
            int[][][][] blockTrace = new int[alphabetSize][][][];
            for (int i = 0; i < alphabetSize; i++)
            {
                // モデルには必ず初期状態が含まれているので、それを除いておく。
                blockLogViterbi[i] = new double[stateNumber - 1][];
                blockTrace[i] = new int[stateNumber - 1][][];
                for (int j = 0; j < stateNumber - 1; j++)
                {
                    blockLogViterbi[i][j] = new double[stateNumber - 1];
                    blockTrace[i][j] = new int[stateNumber - 1][];
                }
            }
            // 初期化
            for (int p = 0; p < alphabetSize; p++)
            {
                for (int i = 0; i < stateNumber - 1; i++)
                {
                    for (int j = 0; j < stateNumber - 1; j++)
                    {
                        blockLogViterbi[p][i][j] = System.Math.Log(emissionProbability[p][i]) + System.Math.Log(transitionProbability[i + 1][j + 1]);
                        blockTrace[p][i][j] = new int[] { i };
                    }
                }
            }
            // 再帰処理
            for (int t = 2; t <= blockLength; t++)
            {
                int subpatternNumber = PowerInt(alphabetSize, t);
                double[][][] previousBlockLogViterbi = blockLogViterbi;
                int[][][][] previousBlockTrace = blockTrace;
                blockLogViterbi = new double[subpatternNumber][][];
                blockTrace = new int[subpatternNumber][][][];
                for (int p = 0; p < subpatternNumber; p++)
                {
                    blockLogViterbi[p] = new double[stateNumber - 1][];
                    blockTrace[p] = new int[stateNumber - 1][][];
                    for (int i = 0; i < stateNumber - 1; i++)
                    {
                        blockLogViterbi[p][i] = new double[stateNumber - 1];
                        blockTrace[p][i] = new int[stateNumber - 1][];
                        for (int j = 0; j < stateNumber - 1; j++)
                        {
                            blockTrace[p][i][j] = new int[t];
                            System.Array.Copy(previousBlockTrace[p / alphabetSize][i][j], blockTrace[p][i][j], t - 1);
                        }
                    }
                }
                for (int p = 0; p < subpatternNumber; p++)
                {
                    int letter = p % alphabetSize;
                    // From state i.
                    for (int i = 0; i < stateNumber - 1; i++)
                    {
                        // To state j.
                        for (int j = 0; j < stateNumber - 1; j++)
                        {
                            double max = double.NegativeInfinity;
                            int argmax = -1;
                            for (int k = 0; k < stateNumber - 1; k++)
                            {
                                double tmp = System.Math.Log(emissionProbability[letter][j]) + System.Math.Log(transitionProbability[k + 1][j + 1]) + previousBlockLogViterbi[p / alphabetSize][i][k];
                                if (tmp > max)
                                {
                                    max = tmp;
                                    argmax = k;
                                }
                            }
                            blockLogViterbi[p][i][j] = max;
                            blockTrace[p][i][j][t - 1] = argmax;
                        }
                    }
                }
            }

            int[] hiddenState = new int[emission.Length + 1];
            // Viterbi variableの計算
            double[] distribution = new double[stateNumber - 1];
            
            int firstCharIndex = System.Array.IndexOf<char>(alphabet, emission[0]);
            for (int i = 0; i < distribution.Length; i++)
            {
                distribution[i] = System.Math.Log(transitionProbability[0][i + 1]) + System.Math.Log(emissionProbability[firstCharIndex][i]);
                hiddenState[0] = -1;
            }
            int fullBlockCount = (emission.Length - 1) / blockLength;
            int[][] trace = new int[fullBlockCount][];
            for (int t = 0; t < fullBlockCount; t++)
            {
                trace[t] = new int[stateNumber - 1];
                // 2文字目から開始。
                int position = t * blockLength + 1;
                double[] previousDistribution = distribution;
                distribution = new double[stateNumber - 1];
                int index = SubstringToIndex(emission, position, blockLength);
                for (int i = 0; i < stateNumber - 1; i++)
                {
                    double max = double.NegativeInfinity;
                    int argmax = -1;
                    for (int j = 0; j < stateNumber - 1; j++)
                    {
                        double tmp = blockLogViterbi[index][j][i] + previousDistribution[j];
                        if (max < tmp)
                        {
                            max = tmp;
                            argmax = j;
                        }
                    }
                    distribution[i] = max;
                    trace[t][i] = argmax;
                }
            }
            

            // 余りの処理
            int[][] traceByOne = new int[(emission.Length - 1) % blockLength][];
            for (int i = 0; i < traceByOne.Length; i++)
            {
                traceByOne[i] = new int[stateNumber - 1];
            }
            for (int position = blockLength * fullBlockCount + 1; position < emission.Length; position++)
            {
                double[] previousDistribution = distribution;
                distribution = new double[stateNumber - 1];
                int charIndex = System.Array.IndexOf<char>(alphabet, emission[position]);
                for (int i = 0; i < stateNumber - 1; i++)
                {
                    double max = double.NegativeInfinity;
                    int argmax = -1;
                    for (int j = 0; j < stateNumber - 1; j++)
                    {
                        double tmp = previousDistribution[i] + System.Math.Log(transitionProbability[i + 1][j + 1]) + System.Math.Log(emissionProbability[charIndex][j]);
                        if (max < tmp)
                        {
                            max = tmp;
                            argmax = j;
                        }
                    }
                    distribution[i] = max;
                    traceByOne[(position - 1) % blockLength][i] = argmax;
                }
            }

            if (debugmode)
            {
                Program.PrintArray<double>(distribution);
            }

            // 終了処理
            double finalMax = double.NegativeInfinity;
            int finalState = -1;
            for (int k = 0; k < stateNumber - 1; k++)
            {
                double tmp = distribution[k];
                if (finalMax < tmp)
                {
                    finalMax = tmp;
                    finalState = k;
                }
            }
            hiddenState[emission.Length] = finalState;
            for (int position = emission.Length - 1; position >= blockLength * fullBlockCount + 1; position--)
            {
                hiddenState[position] = traceByOne[position % blockLength - 1][hiddenState[position + 1]];
            }
            for (int i = trace.Length - 1; i >= 0; i--)
            {
                int index = SubstringToIndex(emission, i * blockLength + 1, blockLength);
                System.Array.Copy(blockTrace[index][trace[i][hiddenState[(i + 1) * blockLength + 1]]][hiddenState[(i + 1) * blockLength + 1]], 0, hiddenState, i * blockLength + 1, blockLength);
            }
            for (int i = 0; i < hiddenState.Length; i++)
            {
                hiddenState[i] += 1;
            }
            return hiddenState;
        }

        /// <summary>
        /// Viterbiアルゴリズムによって、文字列emissionが観測された時の隠れ状態の推定列をint型配列として返す。
        /// </summary>
        /// <param name="emission">観測文字列</param>
        /// <param name="debugMode">[Optional] true時にViterbi変数の対数値、トレースバックテーブルを表示する。既定値はfalse。</param>
        /// <returns>初期状態から始まる、隠れ状態の推定列</returns>
        public int[] Viterbi(string emission, bool debugMode = false)
        {
            // 隠れ状態。これが返り値となる。
            int[] hiddenState = new int[emission.Length + 1];

            // logV means log (ln) of viterbi[time][state]
            double[][] logV = new double[emission.Length + 1][];
            // logVの再帰処理において、logV[t][state] が使用した最大値を与える状態を、trace[t][state]に代入する。
            int[][] trace = new int[emission.Length + 1][];
            for (int t = 0; t <= emission.Length; t++)
            {
                logV[t] = new double[stateNumber];
                trace[t] = new int[stateNumber];
            }

            // t = 0のときの初期化
            logV[0][0] = 0; // ln(1)
            for (int k = 1; k < stateNumber; k++)
            {
                logV[0][k] = double.NegativeInfinity; // ln(0)
            }

            // 再帰処理
            for (int t = 1; t <= emission.Length; t++)
            {
                // t = 0 は初期状態なので、時刻tでの観測文字はemission[t - 1]
                // 観測文字がアルファベット集合の何番目の要素であるかを取得する。
                int currentIndex = System.Array.IndexOf(alphabet, emission[t - 1]);
                for (int l = 0; l < stateNumber; l++)
                {
                    double max = double.NegativeInfinity;
                    int argmax = -1;
                    for (int k = 0; k < stateNumber; k++)
                    {
                        double tmp = logV[t - 1][k] + System.Math.Log(transitionProbability[k][l]);
                        if (max < tmp)
                        {
                            max = tmp;
                            argmax = k;
                        }
                    }
                    if (l == 0)
                    {
                        logV[t][l] = max; // 実際は全て-Infinityとなる。
                    }
                    else
                    {
                        logV[t][l] = System.Math.Log(emissionProbability[currentIndex][l - 1]) + max;
                    }
                    trace[t][l] = argmax;
                }
            }

            // 終了処理
            double finalMax = double.NegativeInfinity;
            int finalState = 0;
            for (int k = 0; k < stateNumber; k++)
            {
                double tmp = logV[emission.Length][k];
                if (finalMax < tmp)
                {
                    finalMax = tmp;
                    finalState = k;
                }
            }

            // 隠れ状態推定列の生成（トレースバック）
            hiddenState[emission.Length] = finalState;
            for (int t = emission.Length - 1; t >= 0; t--)
            {
                hiddenState[t] = trace[t + 1][hiddenState[t + 1]];
            }

            if (debugMode)
            {
                Program.PrintMatrix<double>(logV);
                Program.PrintMatrix<int>(trace);
            }

            return hiddenState;
        }

        /// <summary>
        /// 前向き確率により、観測された文字列が自身のモデルパラメータのもとで観測される確率を返す。
        /// logにtrueが設定されたときは、その確率の自然対数値を返す。
        /// </summary>
        /// <param name="emission">観測文字列</param>
        /// <param name="log">自然対数値を返すかどうかの真偽値</param>
        public double GetForwardProbability(string emission, bool log)
        {
            // スケーリングされた前向き確率のテーブルを取得し、そのうち全確率を求めるのに必要な列を抜き出す。
            double[] scale;
            double[] ret = GetForwardProbability(emission, out scale)[emission.Length];
            double prob = 0;
            for (int i = 0; i < ret.Length; i++)
            {
                prob += ret[i];
            }
            // スケーリングを戻すことで、全確率を計算する。
            // logという変数がtrueのときは、その対数値を計算する。
            if (log)
            {
                prob = System.Math.Log(prob);
                for (int i = 0; i < scale.Length; i++)
                {
                    prob += System.Math.Log(scale[i]);
                }
            }
            else
            {
                for (int i = 0; i < scale.Length; i++)
                {
                    prob *= scale[i];
                }
            }
            return prob;
        }

        /// <summary>
        /// 受け取った観測文字列emissionに対するスケーリングされた前向き確率のテーブルを作成し、
        /// スケーリング　ファクターの配列を引数scaleに格納して返す。
        /// </summary>
        /// <param name="emission">観測文字列</param>
        /// <param name="scale">Scaling factorを格納される配列</param>
        /// <returns>スケーリングされた前向き確率のテーブル</returns>
        double[][] GetForwardProbability(string emission, out double[] scale)
        {
            // 状態k、時刻t（t文字目）におけるスケーリングされた前向き確率をforwardProbability[k][t]に格納する。
            double[][] forwardProbability = new double[emission.Length + 1][];
            for (int i = 0; i <= emission.Length; i++)
            {
                forwardProbability[i] = new double[stateNumber];
            }
            // 配列は自動的に0.0で初期化される。
            scale = new double[emission.Length];
            double[] unscaledProbability = new double[stateNumber];

            // 初期化
            forwardProbability[0][0] = 1;
            for (int i = 1; i < stateNumber; i++)
            {
                forwardProbability[0][i] = 0;
            }

            // 再帰処理
            for (int t = 1; t <= emission.Length; t++)
            {
                // スケーリング前の前向き確率の計算する。
                for (int l = 1; l < stateNumber; l++)
                {
                    double sum = 0;
                    for (int k = 0; k < stateNumber; k++)
                    {
                        sum += forwardProbability[t - 1][k] * transitionProbability[k][l];
                    }
                    unscaledProbability[l] = emissionProbability[System.Array.IndexOf(alphabet, emission[t - 1])][l - 1] * sum;
                }

                // スケーリングファクターの計算
                double currentScale = 0;
                for (int l = 0; l < stateNumber; l++)
                {
                    currentScale += unscaledProbability[l];
                }
                scale[t - 1] = currentScale;

                // スケーリング
                for (int l = 0; l < stateNumber; l++)
                {
                    forwardProbability[t][l] = unscaledProbability[l] / currentScale;
                }
            }
            return forwardProbability;
        }

        /// <summary>
        /// 計算済みの前向き確率の最後列とスケーリングファクターを受け取り、それに対応する全確率の自然対数値を返す。
        /// </summary>
        /// <param name="lastForwardProbability">計算済みの前向き確率の最後列(t = T)</param>
        /// <param name="scale">スケーリングファクター</param>
        /// <returns>全確率の自然対数値</returns>
        double UnscaleForwardProbabilityLog(double[] lastForwardProbability, double[] scale)
        {
            double prob = 0;
            for (int i = 0; i < lastForwardProbability.Length; i++)
            {
                prob += lastForwardProbability[i];
            }
            prob = System.Math.Log(prob);
            for (int i = 0; i < scale.Length; i++)
            {
                prob += System.Math.Log(scale[i]);
            }
            return prob;
        }

        /// <summary>
        /// 後向き確率により、観測された文字列が自身のモデルパラメータのもとで観測される確率を返す。
        /// logにtrueが設定されたときは、その確率の自然対数値を返す。
        /// </summary>
        /// <param name="emission">観測文字列</param>
        /// <param name="log">自然対数値を返すかどうかの真偽値</param>
        public double GetBackwardProbability(string emission, bool log)
        {
            double[] scale;
            double[] ret = GetBackwardProbability(emission, out scale)[1];
            double prob = 0;
            for (int i = 1; i < ret.Length; i++)
            {
                prob += transitionProbability[0][i] * emissionProbability[System.Array.IndexOf(alphabet, emission[0])][i - 1] * ret[i];
            }
            // スケーリングを戻す。
            if (log)
            {
                prob = System.Math.Log(prob);
                for (int i = 0; i < scale.Length - 1; i++)
                {
                    prob += System.Math.Log(scale[i]);
                }
            }
            else
            {
                for (int i = 0; i < scale.Length - 1; i++)
                {
                    prob *= scale[i];
                }
            }
            return prob;
        }

        /// <summary>
        /// 受け取った観測文字列emissionに対するスケーリングされた後向き確率のテーブルを作成し、
        /// スケーリング　ファクターの配列を引数scaleに格納して返す。
        /// </summary>
        /// <param name="emission">観測文字列</param>
        /// <param name="scale">Scaling factorを格納される配列</param>
        /// <returns>スケーリングされた後向き確率のテーブル</returns>
        double[][] GetBackwardProbability(string emission, out double[] scale)
        {
            double[][] backwardProbability = new double[emission.Length + 1][];
            for (int i = 0; i <= emission.Length; i++)
            {
                backwardProbability[i] = new double[stateNumber];
            }
            // 配列は自動的に0.0で初期化される。
            scale = new double[emission.Length];
            double[] unscaledProbability = new double[stateNumber];

            // 初期化
            backwardProbability[emission.Length][0] = 0;
            for (int i = 1; i < stateNumber; i++)
            {
                backwardProbability[emission.Length][i] = 1;
            }

            // 再帰処理
            for (int t = emission.Length - 1; t >= 1; t--)
            {
                // スケーリング前の後向き確率の計算。
                for (int k = 1; k < stateNumber; k++)
                {
                    double sum = 0;
                    for (int l = 1; l < stateNumber; l++)
                    {
                        sum += transitionProbability[k][l] * backwardProbability[t + 1][l] * emissionProbability[System.Array.IndexOf(alphabet, emission[t])][l - 1];
                    }
                    unscaledProbability[k] = sum;
                }
                // スケーリングファクターの計算。
                double currentScale = 0;
                for (int k = 0; k < stateNumber; k++)
                {
                    currentScale += unscaledProbability[k];
                }
                scale[t - 1] = currentScale;
                // スケーリング
                for (int k = 0; k < stateNumber; k++)
                {
                    backwardProbability[t][k] = unscaledProbability[k] / currentScale;
                }
            }
            return backwardProbability;
        }

        /// <summary>
        /// Baum-Welchアルゴリズムにより、与えられた観測文字列の配列emissionsを学習データに用い、
        /// 自身の遷移確率行列と出力確率行列を、現在の値を初期値として、モデルの尤度が局所最大値（の近似値）を取るときの値に更新する。
        /// </summary>
        /// <param name="emissions">観測文字列の配列</param>
        /// <param name="delta">終了条件となる対数尤度の変化量の境界値</param>
        /// <param name="verboseMode">[Optional] 再帰処理の各ステップでの対数尤度の表示設定。既定値はfalse</param>
        public void BaumWelch(string[] emissions, double delta, bool verboseMode = false)
        {
            // 遷移や出力の頻度を数える（期待値計算）ための2次元配列を用意する。
            // transitionFrequency["from state"]["to state"]
            double[][] transitionFrequency = new double[stateNumber][];
            for (int i = 0; i < stateNumber; i++)
            {
                transitionFrequency[i] = new double[stateNumber];
            }
            // emissionFrequency["state excluding initial state"]["character index"]
            double[][] emissionFrequency = new double[stateNumber - 1][];
            for (int i = 0; i < stateNumber - 1; i++)
            {
                emissionFrequency[i] = new double[alphabetSize];                
            }

            bool notFirst = false;
            double previousLogLikelihood = double.NegativeInfinity;
            while (true)
            {
                // 遷移や出力の期待値を計算する。同時に期待値計算に利用したパラメータにおける対数尤度を返す。
                double currentLogLikelihood = CalculateExpectation(emissions, transitionFrequency, emissionFrequency);
                // 対数尤度の変化量を求める。
                double diff = currentLogLikelihood - previousLogLikelihood;
                // verboseModeがtrueなら、対数尤度の変化を表示する。
                if (verboseMode)
                {
                    if (notFirst)
                    {
                        // 下記は初回のみ実行されない。
                        Console.WriteLine("               > {0:F15}", diff);
                    }
                    Console.WriteLine("{0:F15}", currentLogLikelihood);
                }
                notFirst = true;
                // 対数尤度の変化量diffがdeltaより小さいとき、アルゴリズムを終了する。
                if (diff > delta)
                {
                    // 期待値をもとにモデルのパラメータを更新する。
                    UpdateParameters(transitionFrequency, emissionFrequency);
                    // 対数尤度の更新。
                    previousLogLikelihood = currentLogLikelihood;
                }
                else
                {
                    break;
                }
            }
            return;
        }

        /// <summary>
        /// 受け取った学習データemissionsに対して、現在のモデルパラメータを用いて状態遷移と出力記号の期待値を計算し、
        /// それぞれtransitionFrequencyとemissionFrequencyに格納する。
        /// 返り値として、現在のモデルパラメータによる対数尤度を返す。
        /// </summary>
        /// <param name="emissions">学習データ（観測文字列の配列）</param>
        /// <param name="transitionFrequency">状態遷移の期待値を格納する2次元配列</param>
        /// <param name="emissionFrequency">出力記号の期待値を格納する2次元配列</param>
        /// <returns>計算に用いたモデルパラメータによる対数尤度</returns>
        private double CalculateExpectation(string[] emissions, double[][] transitionFrequency, double[][] emissionFrequency)
        {
            double logLikelihood = 0.0;
            // 学習データを1つずつ取り出して期待値を計算し、配列の対応する箇所に加算していく。
            for (int j = 0; j < emissions.Length; j++)
            {
                // 現在の学習配列に対する前向き確率、後向き確率とそれぞれのスケーリングファクターを得る。
                double[] forwardScale;
                double[][] forwardProbability = GetForwardProbability(emissions[j], out forwardScale);
                double[] backwardScale;
                double[][] backwardProbability = GetBackwardProbability(emissions[j], out backwardScale);
                // 状態遷移の期待値を計算する。
                for (int k = 0; k < stateNumber; k++)
                {
                    for (int l = 1; l < stateNumber; l++)
                    {
                        double sum = 0.0;
                        for (int t = 0; t < emissions[j].Length; t++)
                        {
                            double product = 1.0;
                            // スケーリングの解除のためにスケーリングファクターをまとめて計算する。
                            for (int i = t; i < backwardScale.Length - 1; i++)
                            {
                                product = product * backwardScale[i] / forwardScale[i];
                            }
                            product /= forwardScale[forwardScale.Length - 1];
                            sum += forwardProbability[t][k] * transitionProbability[k][l] * emissionProbability[System.Array.IndexOf(alphabet, emissions[j][t])][l - 1] * backwardProbability[t + 1][l] * product;
                        }
                        // 配列は作成時に初期値（0.0）で初期化されている。
                        transitionFrequency[k][l] += sum;
                    }
                }
                // 出力記号の期待値を計算する。
                for (int k = 0; k < stateNumber - 1; k++)
                {
                    for (int b = 0; b < alphabetSize; b++)
                    {
                        double sum = 0.0;
                        for (int t = 1; t <= emissions[j].Length; t++)
                        {
                            if (b == System.Array.IndexOf(alphabet, emissions[j][t - 1]))
                            {
                                double product = 1.0;
                                // スケーリング解除のためにスケーリングファクターをまとめて計算する。
                                for (int i = t; i < backwardScale.Length - 1; i++)
                                {
                                    product = product * backwardScale[i] / forwardScale[i];
                                }
                                product = product * backwardScale[t - 1] / forwardScale[forwardScale.Length - 1];
                                sum += forwardProbability[t][k + 1] * backwardProbability[t][k + 1] * product;
                            }
                        }
                        emissionFrequency[k][b] += sum;
                    }
                }
                // 現在の学習配列に対する対数尤度（全確率の自然対数値）を求め、全体の対数尤度に加算する。
                logLikelihood += UnscaleForwardProbabilityLog(forwardProbability[forwardProbability.Length - 1], forwardScale);
            }
            return logLikelihood;
        }

        /// <summary>
        /// 計算された期待値をもとにモデルパラメータを更新する。
        /// 更に期待値を格納していた2次元配列の各要素を0.0で初期化する。
        /// </summary>
        /// <param name="transitionFrequency">状態遷移の期待値を格納した2次元配列</param>
        /// <param name="emissionFrequency">出力記号の期待値を格納した2次元配列</param>
        private void UpdateParameters(double[][] transitionFrequency, double[][] emissionFrequency)
        {
            // 期待値から状態遷移確率を求める。
            for (int k = 0; k < stateNumber; k++)
            {
                double sum = 0.0;
                for (int l = 0; l < stateNumber; l++)
                {
                    sum += transitionFrequency[k][l];
                }
                for (int l = 0; l < stateNumber; l++)
                {
                    // 状態遷移確率の更新。
                    transitionProbability[k][l] = transitionFrequency[k][l] / sum;
                    // 次の計算のために初期化。
                    transitionFrequency[k][l] = 0.0;
                }
            }
            // 期待値から出力確率を求める。
            for (int k = 0; k < stateNumber - 1; k++)
            {
                double sum = 0.0;
                for (int b = 0; b < alphabetSize; b++)
                {
                    sum += emissionFrequency[k][b];
                }
                for (int b = 0; b < alphabetSize; b++)
                {
                    // 出力確率の更新。
                    emissionProbability[b][k] = emissionFrequency[k][b] / sum;
                    // 次の計算のために初期化する。
                    emissionFrequency[k][b] = 0.0;
                }
            }
        }

    }
}
