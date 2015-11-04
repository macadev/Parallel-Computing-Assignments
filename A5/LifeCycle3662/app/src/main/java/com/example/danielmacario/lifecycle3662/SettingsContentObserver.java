package com.example.danielmacario.lifecycle3662;

import android.content.Context;
import android.database.ContentObserver;
import android.media.AudioManager;
import android.net.Uri;
import android.os.Handler;
import android.util.Log;
import android.widget.Toast;

public class SettingsContentObserver extends ContentObserver {
    private AudioManager audioManager;

    private String TAG = "Volume Content Observer";
    private Context context;

    public SettingsContentObserver(Context context, Handler handler) {
        super(handler);
        this.context = context;
        audioManager = (AudioManager) context.getSystemService(Context.AUDIO_SERVICE);
    }

    @Override
    public boolean deliverSelfNotifications() {
        return false;
    }

    @Override
    public void onChange(boolean selfChange) {
        int currentVolume = audioManager.getStreamVolume(AudioManager.STREAM_SYSTEM);
        int duration = Toast.LENGTH_SHORT;

        Toast toast = Toast.makeText(context, "Volume Changed! " + currentVolume, duration);
        toast.show();
    }

}