using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace BRSSetup
{
    public partial class Form1 : Form
    {
        public string selected="";
        public Form1()
        {
            InitializeComponent();
        }
        private void button1_Click(object sender, EventArgs e)
        {
            string homePath=System.IO.Path.GetDirectoryName(Application.ExecutablePath).ToString();
            char[] name = selected.ToCharArray();
            string fileName = new String(name, 0, 3)+ name[4].ToString()+ name[name.Length-1].ToString();

            string fromPath = homePath + "/configs/" + fileName +".config";
            string toPath = homePath + "/test.config";
            File.Copy(fromPath, toPath);
        }

        private void listBox2_SelectedIndexChanged(object sender, EventArgs e)
        {
            listBox1.ClearSelected();
            try
            {
                selected = listBox2.SelectedItem.ToString();
            }
            catch(Exception ex)
            {
            }
        }

        private void listBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            listBox2.ClearSelected();
            try
            {
                selected = listBox1.SelectedItem.ToString();
            }
            catch (Exception ex)
            {
            }
        }
        
    }
}
