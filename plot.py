import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def extract_data(runname):
    df = pd.read_csv(
        runname + "/analysis.txt",
        sep=r",\s*",
        engine="python",
        header=None,
        names=["time", "Etot", "Rvir"]
    )

    df = df.replace(
        {r"^[^=]*=\s*": ""},
        regex=True
    ).astype(float)

    data = df.to_numpy()
    print(data.shape)

    return data

def main():
    runs =["1","2","4","8","16",]

    datas = [extract_data(run) for run in runs]


    # Figureの初期化
    fig = plt.figure(figsize=(12, 8)) #...1

    # Figure内にAxesを追加()
    ax = fig.add_subplot(211) #...2

    ax.set_xlabel("time")
    ax.set_ylabel("etot")

    ax2 = fig.add_subplot(212) #...2
    
    ax2.set_xlabel("time")
    ax2.set_ylabel("rvir")
    for i in range(len(runs)):    
        ax.plot(datas[i][:,0], datas[i][:,1], label = runs[i]) #...3
        ax2.plot(datas[i][:,0], datas[i][:,2], label = runs[i]) #...3
    # 凡例の表示
    plt.legend()

    # プロット表示(設定の反映)
    plt.savefig("out.png")

    errors=[]
    resolution=[]
    for i in range(len(datas)):
        error = datas[i][-1,1] - datas[i][0,1]
        errors.append(error)
        resolution.append(1/float(runs[i]))

    fig1= plt.figure(figsize=(12,8))
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel("Δt/default Δt")
    ax1.set_ylabel("total energy error")

    ax1.scatter(resolution, errors)

    
    res2=np.polyfit(resolution, errors, 2)
    print(res2)
    x_plot = np.linspace(0.,1.,100)
    y2 = np.poly1d(res2)(x_plot)
    # プロット表示(設定の反映)
    ax1.plot(x_plot,y2)
    plt.savefig("error.png")




if __name__ == '__main__':
    main()

