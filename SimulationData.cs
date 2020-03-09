using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;

namespace Simulation.FEM
{
    class SimulationData
    {

        public class Inner
        {
            public double conductivity;
            public double specificHeat;
            public double ro;
            public int alfa;
            public double ambientTemperature; //otoczenia
            public double simulationTime;
            public double dTau;
            public double initialtemp;
            public int nodesHeight;
            public int nodesLength;
            public double realHeight;
            public double realLength;
        }

        private static SimulationData _instance;
        private Inner inner;

        public static SimulationData Instance
        {
            get
            {
                return _instance;
            }
        }

        public double Conductivity { get => inner.conductivity; }
        public double SpecificHeat { get => inner.specificHeat;  }
        public double Ro { get => inner.ro;  }
        public int Alfa { get => inner.alfa; }
        public double AmbientTemperature { get => inner.ambientTemperature; }
        public double Time { get => inner.simulationTime; }
        public double DTau { get => inner.dTau;  }
        public double Initialtemp { get => inner.initialtemp; }
        public int NodesHeight { get => inner.nodesHeight; }
        public int NodesLength { get => inner.nodesLength; }
        public double RealHeight { get => inner.realHeight; }
        public double RealLength { get => inner.realLength; }


        public static void LoadFromFile(string path)
        {
            string file = File.ReadAllText(path);
            Inner deserializedInner = JsonConvert.DeserializeObject<Inner>(file);
            _instance = new SimulationData();
            _instance.inner = deserializedInner;
        }

       
    }
}
