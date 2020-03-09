using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace Simulation
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
          
            try
            {
                FEM.SimulationData.LoadFromFile("./SimulationData2.json");
            }
            catch (IOException) 
            {
                Console.WriteLine("File not found");
                return; 
            }

            FEM.GlobalData data = new FEM.GlobalData { NHeight = FEM.SimulationData.Instance.NodesHeight,
                                                       NLength = FEM.SimulationData.Instance.NodesLength, 
                                                       GlobalHeightGrid = FEM.SimulationData.Instance.RealHeight,
                                                       GlobalLengthGrid = FEM.SimulationData.Instance.RealLength };
            FEM.Grid grid = new FEM.Grid(data);
            grid.ConstructGrid();
            
            
            
        }
    }
}