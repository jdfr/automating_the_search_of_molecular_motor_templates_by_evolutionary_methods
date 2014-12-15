package jd;


import ch.ethz.ssh2.InteractiveCallback;

/**
 * This class is a very simple fake to use the keyboard-interactive
 * authentication method without prompting the password to the user, 
 * but doing it programmatically. Of course, it cannot handle dynamic
 * challenges, only predefined ones 
 * 
 * @author Josedavid
 *
 */
public class PreloadedInteractive implements InteractiveCallback {

	
	    public int count;
        public int r;
		public boolean showThings;
		public String[] responses;

		public PreloadedInteractive(String[] responses, boolean showThings)
		{
			count = 0;
			r = 0;
			this.responses  = responses;
			this.showThings = showThings;
		}
		public PreloadedInteractive()
		{
			count = 0;
			r = 0;
			responses  = null;
			showThings = true;
		}
		
		public String[] getResponses() {
			String[] results = new String[responses.length];
			for (int i = 0; i < results.length; i++) {
				results[i] = responses[i];
			}
			  return results;
			}
		public void setResponses(String[] responses) {
			r = 0;
			this.responses  = responses;
			}
		public int getCursor() {
			  return r;
			}
		public void setCursor(int r) {
			this.r = r;
			}
		public boolean getShowThings() {
			  return showThings;
			}
		public void setShowThings(boolean showThings) {
			this.showThings = showThings;
			}
		public void resetCount() {
			count = 0;
		}

		public String[] replyToChallenge(String name, String instruction, int numPrompts, String[] prompt,
				boolean[] echo) 
		{
			count++;
			if (showThings) {
			  System.out.println("\nntimes="+count+"name=<<"+name+">>\ninstruction=<<"+instruction+">>\nnumPrompts="+numPrompts);
			  for (int i = 0; i < prompt.length; i++) {
				System.out.println("prompt "+i+": echo="+echo[i]+", <<"+prompt[i]+">>");
			  }
			  System.out.println("");
			}
			  
			String[] result = new String[numPrompts];

			if (responses==null) {
				for (int i = 0; i < result.length; i++) {
					result[i] = "";
				}
			} else if (responses.length==0) {
				for (int i = 0; i < result.length; i++) {
					result[i] = "";
				}
			} else for (int i = 0; i < numPrompts; i++)
			{
				result[i] = responses[r++];
				if (r>=responses.length) {
					r = 0;
				}
			}

			return result;
		}

}
