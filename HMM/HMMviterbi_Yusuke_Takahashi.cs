using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

namespace HMM
{
    class Program
    {
        /// <summary>
        /// 1つ目のファイルパスからHMMのパラメータを読み取り、HMMを作成し、2つ目のfasta形式のファイルパスから観測文字列を読み取る。
        /// このHMMに観測文字列を与えることで、隠れ状態の推定などを行う。
        /// </summary>
        /// <param name="args">コマンドライン引数</param>
        static void Main(string[] args)
        {
            //Console.WriteLine(double.NegativeInfinity);

            HMM model = new HMM(0, null, 1);

            // ファイル読み出しのストリームを用意。
            System.IO.StreamReader stream = null;
            // try-catch-finally文
            // 例外発生時にcatch文が実行され、例外の有無に関わらずfinally文が実行される。
            try
            {
                stream = new System.IO.StreamReader(args[0]);
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
                Console.WriteLine("param:{0}", e.Message);
                return;
            }
            finally
            {
                // ファイル読み込みの成功に関わらず実行する処理
                if (stream != null)
                {
                    stream.Close();
                }
            }

            // 観測文字列をファイルから読み込む。
            // string型より文字列操作において高速なStringBuilderクラスを利用する。
            const int defaultCapacity = 100000;
            System.Text.StringBuilder seq = new System.Text.StringBuilder(defaultCapacity);
            try
            {
                stream = new System.IO.StreamReader(args[1]);
                // FASTA形式ファイルの最初の1行を読み飛ばす。
                stream.ReadLine();
                while (!stream.EndOfStream)
                {
                    seq.Append(stream.ReadLine());
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("sequence:{0}", e.Message);
                return;
            }
            finally
            {
                if (stream != null)
                {
                    stream.Close();
                }
            }

            // Viterbiアルゴリズムを読み取った文字列に対して行う。
            PrintArray(model.Viterbi(seq.ToString()));

            /*
            Console.WriteLine(model.GetForwardProbability(seq.ToString(), true).ToString("G17"));
            Console.WriteLine(model.GetBackwardProbability(seq.ToString(), true).ToString("G17"));
            Console.ReadKey(true);
            */
            return;
        }

        /// <summary>
        /// 配列の要素を半角スペース区切りで表示する。
        /// </summary>
        static void PrintArray(int[] array)
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
                    Console.Write("{0} ", matrix[i][j]);
                }
                Console.WriteLine("");
            }
            return;
        }
    }

    /// <summary>
    /// HMMのパラメータを保有し、各種メソッドを提供するクラス
    /// </summary>
    class HMM
    {
        // コンストラクタ以降に変更できないメンバ変数を用意する。
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

        /// <summary>
        /// Viterbiアルゴリズムによって、文字列emissionが観測された時の隠れ状態の推定列をint型配列として返す。
        /// </summary>
        /// <param name="emission">観測文字列</param>
        /// <returns>初期状態から始まる、隠れ状態の推定列</returns>
        public int[] Viterbi(string emission)
        {
            // デバッグ情報の表示
            const bool debugMode = false;
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
        // Viterbiはここまで。

        public double GetForwardProbability(string emission, bool log)
        {
            double[] scale;
            double[] ret = GetForwardProbability(emission, out scale)[emission.Length];
            double prob = 0;
            for (int i = 0; i < ret.Length; i++)
            {
                prob += ret[i];
            }
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

        double[][] GetForwardProbability(string emission, out double[] scale)
        {
            double[][] forwardProbability = new double[emission.Length + 1][];
            for (int i = 0; i <= emission.Length; i++)
            {
                forwardProbability[i] = new double[stateNumber];
            }
            // 配列は自動的に0.0で初期化される。
            scale = new double[emission.Length];
            double[] unscaledProbability = new double[stateNumber];

            forwardProbability[0][0] = 1;
            for (int i = 1; i < stateNumber; i++)
            {
                forwardProbability[0][i] = 0;
            }

            for (int t = 1; t <= emission.Length; t++)
            {
                for (int l = 1; l < stateNumber; l++)
                {
                    double sum = 0;
                    for (int k = 0; k < stateNumber; k++)
                    {
                        sum += forwardProbability[t - 1][k] * transitionProbability[k][l];
                    }
                    unscaledProbability[l] = emissionProbability[System.Array.IndexOf(alphabet, emission[t - 1])][l - 1] * sum;
                }
                double currentScale = 0;
                for (int l = 0; l < stateNumber; l++)
                {
                    currentScale += unscaledProbability[l];
                }
                /*
                // スケール無し。 
                currentScale = 1.0;
                */
                scale[t - 1] = currentScale;
                for (int l = 0; l < stateNumber; l++)
                {
                    forwardProbability[t][l] = unscaledProbability[l] / currentScale;
                }
            }
            return forwardProbability;
        }

        public double GetBackwardProbability(string emission, bool log)
        {
            double[] scale;
            double[] ret = GetBackwardProbability(emission, out scale)[1];
            double prob = 0;
            for (int i = 1; i < ret.Length; i++)
            {
                prob += transitionProbability[0][i] * emissionProbability[System.Array.IndexOf(alphabet, emission[0])][i - 1] * ret[i];
            }
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

            backwardProbability[emission.Length][0] = 0;
            for (int i = 1; i < stateNumber; i++)
            {
                backwardProbability[emission.Length][i] = 1;
            }

            for (int t = emission.Length - 1; t >= 1; t--)
            {
                for (int k = 1; k < stateNumber; k++)
                {
                    double sum = 0;
                    for (int l = 1; l < stateNumber; l++)
                    {
                        sum += transitionProbability[k][l] * backwardProbability[t + 1][l] * emissionProbability[System.Array.IndexOf(alphabet, emission[t])][l - 1];
                    }
                    unscaledProbability[k] = sum;
                }
                double currentScale = 0;
                for (int k = 0; k < stateNumber; k++)
                {
                    currentScale += unscaledProbability[k];
                }
                /*
                // スケール無し。 
                currentScale = 1.0;
                */
                scale[t - 1] = currentScale;
                for (int k = 0; k < stateNumber; k++)
                {
                    backwardProbability[t][k] = unscaledProbability[k] / currentScale;
                }
            }
            return backwardProbability;
        }

        /// <summary>
        /// Baum-Welchアルゴリズムにより、与えられた観測文字列の配列emissionsを用いて、
        /// 自身の遷移確率行列と出力確率行列を、モデルの尤度が局所最大値（の近似値）を取るときの値に書き換える。
        /// </summary>
        /// <param name="emissions">観測文字列の配列</param>
        public void BaumWelch(string[] emissions)
        {
            double[][] transitionSum = new double[stateNumber][];
            double[][] emissionSum = new double[stateNumber][];
            for (int i = 0; i < stateNumber; i++)
            {
                transitionSum[i] = new double[stateNumber];
                emissionSum[i] = new double[alphabetSize];
            }
            for (int k = 0; k < stateNumber; k++)
            {
                for (int l = 0; l < stateNumber; l++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < emissions.Length; j++)
                    {
                        double[] forwardScale;
                        double[][] forwardProbability = GetForwardProbability(emissions[j], out forwardScale);
                        double[] backwardScale;
                        double[][] backwardProbability = GetBackwardProbability(emissions[j], out backwardScale);
                        for (int t = 0; t < emissions[j].Length; t++)
                        {
                            for (int i = 0; i < forwardScale.Length; i++)
                            {
                                
                            }
                        }
                    }
                }
            }
            return;
        }

    }
}
