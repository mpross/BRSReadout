using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BRSReadout
{
    class EmailError
    {
        public static void emailAlert(Exception ex)
        {
            StringBuilder words = new StringBuilder();

            String toEmail = "mpross2@uw.edu";
            String fromEmail = "beamrotation@gmail.com";
            String user = "beamrotation";
            String pass = "windyTilt";

            words.Append("Message: ");
            words.Append(ex.Message);
            words.Append("\n");
            words.Append("Source: ");
            words.Append(ex.Source);
            words.Append("\n");
            words.Append("Stack Trace: ");
            words.Append(ex.StackTrace);
            words.Append("\n");
            words.Append("Target: ");
            words.Append(ex.TargetSite);
            words.Append("\n");

            System.Net.Mail.MailMessage message = new System.Net.Mail.MailMessage();
            message.To.Add(toEmail);
            message.Subject = "LLO Input-Y BRS Error Alert";
            message.From = new System.Net.Mail.MailAddress(fromEmail);
            message.Body = words.ToString();
            System.Net.Mail.SmtpClient smtp = new System.Net.Mail.SmtpClient("smtp.gmail.com");
            smtp.Port = 25;
            smtp.EnableSsl = true;
            smtp.Credentials = new System.Net.NetworkCredential(user, pass);
            //smtp.Send(message);

        }
    }
}
