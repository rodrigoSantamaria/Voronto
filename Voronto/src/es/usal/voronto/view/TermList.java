package es.usal.voronto.view;

import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JButton;
import javax.swing.ListModel;
import javax.swing.JScrollPane;

import es.usal.voronto.control.Voronto;

import java.awt.Color;

public class TermList extends JFrame{
	public Voronto parent=null;
	private Introduction intro=null;
	private JList list;
	private JButton btnOk;
	private JList jlist;
	ArrayList<String> names=null;

	
	public TermList(ArrayList<String> names, Introduction intro, Voronto parent) {
		super();
		this.names=names;
		this.parent=parent;
		this.intro=intro;
		initGUI();
	}
	
	public void initGUI()
		{
		getContentPane().setBackground(Color.WHITE);
		this.setSize(350, 350);
		this.parent=parent;
		setTitle("Select root term");
		
		getContentPane().setLayout(null);
		
		btnOk = new JButton("Ok");
		btnOk.setBackground(Color.WHITE);
		btnOk.setBounds(40, 293, 117, 29);
		getContentPane().add(btnOk);
		this.getRootPane().setDefaultButton(btnOk);
		btnOk.requestFocus();
		btnOk.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) 
				{
				parent.root=(String)list.getSelectedValue();
				intro.jLabel9.setText("Root (click to select): "+parent.root);
				setVisible(false);
				}
			});
		
		JButton btnCancel = new JButton("Cancel");
		btnCancel.setBounds(197, 293, 117, 29);
		getContentPane().add(btnCancel);
		btnCancel.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) 
				{
				setVisible(false);
				}
			});
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(6, 6, 338, 282);
		getContentPane().add(scrollPane);
		
		list = new JList();
		ListModel jListModel = 
				new DefaultComboBoxModel(names.toArray(new String[0]));
		list.setModel(jListModel);
		list.setSelectedIndex(0);
		list.setAutoscrolls(true);
		scrollPane.setViewportView(list);
		}
}
