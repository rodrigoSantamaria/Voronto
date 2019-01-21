package es.usal.voronto.view;

import java.awt.event.KeyEvent;

import processing.core.PApplet;
import processing.core.PFont;

public class PSearchFrame extends PApplet{
	private static final long serialVersionUID = 3670672800541188702L;
	private VoronoiVisualization vv;
	public String searchText="";
	private PFont font;
	private int fontSize;
	public PSearchFrame(VoronoiVisualization v)
		{
		super();
		
		vv=v;
		width=200;
		height=50;
		font=vv.font;
		fontSize=vv.fontSize;
		setSize(width, height);
		}

	public void setup()
		{
		smooth();
		noLoop();
		}
	
	public void draw()
		{
		fill(255,255,255);
		stroke(0);
		strokeWeight(1);
		rect(0,0,width,height);
		
		fill(154);
		textFont(font, fontSize);
		textAlign(CENTER, TOP);
		text("Search for:",(int)(width*0.5),5);
		text(searchText,(int)(width*0.5),25);
		}
	
	public void keyReleased()
		{
		switch(keyCode)
			{
			case KeyEvent.VK_ESCAPE:
				vv.search(null);
				break;
			//TODO: escape
			case 16://others (to ignore)
				return;
			case 10://enter
				vv.search(searchText);
				return;
				
			case 8://supr
				if(searchText.length()>0)	
					searchText=searchText.substring(0, searchText.length()-1);
				redraw();
				return;
			default:
				searchText=searchText+""+key;
				redraw();
				return;
			}
		}
	public void exit()
		{
		System.out.println("Exiting...");
		}

}
